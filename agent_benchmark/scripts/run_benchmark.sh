#!/usr/bin/env bash
# Agent Benchmark: 完全隔离运行脚本
# 用法: run_benchmark.sh <agent>
#   agent: piggpa | claudecode | openclaw
#
# 隔离策略:
# - claudecode/openclaw: 使用 bwrap 沙箱，仅挂载 data + 自身工作区 + 工具
# - piggpa: 使用自身内部沙箱，cwd = 自身工作区，无限时
# - 三方均清除持久状态（无记忆、无 skills、无配置）
#
# NOTE: This copy has been redacted for public release. Real API keys / gateway
# tokens / absolute paths used at run time have been replaced with placeholders.
set -u

AGENT="$1"
ROOT="/workspace/pigbole/benchmark/agent_benchmark"
TASK_ID="agent_benchmark"
TASK_DIR="${ROOT}/task"
DATA_DIR="${ROOT}/data"
LOG_DIR="${ROOT}/logs"

# 每个 agent 的工作区
case "${AGENT}" in
  piggpa)     AGENT_DIR="${ROOT}/piggpa" ;;
  claudecode) AGENT_DIR="${ROOT}/claudecode" ;;
  openclaw)   AGENT_DIR="${ROOT}/openclaw" ;;
  *) echo "Unknown agent: ${AGENT}" >&2; exit 1 ;;
esac

mkdir -p "${AGENT_DIR}" "${LOG_DIR}"

# 读取 prompt（统一，无 agent 特定替换）
PROMPT=$(cat "${TASK_DIR}/prompt.txt")
echo "${PROMPT}" > "${AGENT_DIR}/prompt.txt"

START_TS=$(date +%s)
START_ISO=$(date -Iseconds)

SESSION_LOG="${AGENT_DIR}/session_log.txt"
CONTEXT_JSON="${AGENT_DIR}/session_context.json"
METRICS_JSON="${AGENT_DIR}/metrics.json"
SESSION_ID=""
EXIT_CODE=0
IN_TOK=0; OUT_TOK=0; CACHE_R=0; CACHE_W=0; REAS_TOK=0

# ============================================================
# 清除持久状态（无记忆）
# ============================================================
clear_state() {
  local ts=$(date +%s)
  local trash="/tmp/ai_trash_${ts}"
  mkdir -p "${trash}"

  # Claude Code 持久状态
  if [ -d "${HOME}/.claude/projects" ]; then
    mv "${HOME}/.claude/projects" "${trash}/claude_projects" 2>/dev/null || true
    mkdir -p "${HOME}/.claude/projects"
  fi
  if [ -d "${HOME}/.claude/sessions" ]; then
    mv "${HOME}/.claude/sessions"/* "${trash}/" 2>/dev/null || true
  fi

  # OpenClaw 持久状态
  if [ -d "${HOME}/.openclaw/state" ]; then
    mv "${HOME}/.openclaw/state"/* "${trash}/" 2>/dev/null || true
  fi
  if [ -d "${HOME}/.openclaw/cache" ]; then
    mv "${HOME}/.openclaw/cache"/* "${trash}/" 2>/dev/null || true
  fi

  echo "[state] cleared persistent state for ${AGENT} (moved to ${trash})"
}

# ============================================================
# bwrap 沙箱公共挂载
# ============================================================
bwrap_common_args() {
  echo \
    --ro-bind /usr /usr \
    --ro-bind /lib /lib \
    --ro-bind /lib64 /lib64 \
    --ro-bind /bin /bin \
    --ro-bind /sbin /sbin \
    --ro-bind /etc /etc \
    --dev /dev \
    --proc /proc \
    --ro-bind /sys /sys \
    --ro-bind /workspace/miniconda3 /workspace/miniconda3 \
    --ro-bind /workspace/app /workspace/app \
    --ro-bind /workspace/software /workspace/software \
    --ro-bind "${DATA_DIR}" "${DATA_DIR}" \
    --bind "${AGENT_DIR}" "${AGENT_DIR}" \
    --bind /tmp /tmp
}

# ============================================================
# 运行各 Agent
# ============================================================
case "${AGENT}" in
  piggpa)
    # PigGPA: 无限时，使用自身内部沙箱，每次 chat -q 创建新 session（无记忆）
    PIGGPA_RAW="${LOG_DIR}/piggpa_raw_$$.log"
    cd "${AGENT_DIR}"
    # 注意: 不设 --max-turns（无限时）
    piggpa chat -q "${PROMPT}" > "${PIGGPA_RAW}" 2>&1
    EXIT_CODE=$?
    cp "${PIGGPA_RAW}" "${SESSION_LOG}"

    # 提取 session_id
    SESSION_ID=$(grep -E '^session_id:' "${PIGGPA_RAW}" | tail -1 | awk '{print $2}')
    if [ -z "${SESSION_ID}" ]; then
      # 尝试从 piggpa session 文件找最新的
      SESSION_ID=$(ls -t /workspace/pigbole/.piggpa/sessions/session_*.json 2>/dev/null | head -1 | sed 's|.*/session_||; s|\.json||')
    fi

    # 拷贝 session JSON
    if [ -n "${SESSION_ID}" ] && [ -f "/workspace/pigbole/.piggpa/sessions/session_${SESSION_ID}.json" ]; then
      cp "/workspace/pigbole/.piggpa/sessions/session_${SESSION_ID}.json" "${CONTEXT_JSON}"
    else
      echo '{"error": "no session json"}' > "${CONTEXT_JSON}"
    fi

    # Token from sqlite
    if [ -n "${SESSION_ID}" ]; then
      TOKEN_ROW=$(sqlite3 /workspace/pigbole/.piggpa/state.db "SELECT COALESCE(input_tokens,0), COALESCE(output_tokens,0), COALESCE(cache_read_tokens,0), COALESCE(cache_write_tokens,0), COALESCE(reasoning_tokens,0) FROM sessions WHERE id='${SESSION_ID}';" 2>/dev/null || echo "0|0|0|0|0")
      IFS='|' read -r IN_TOK OUT_TOK CACHE_R CACHE_W REAS_TOK <<< "${TOKEN_ROW}"
    fi
    ;;

  claudecode)
    # Claude Code: bwrap 沙箱，无持久状态
    CC_RAW="${LOG_DIR}/cc_raw_$$.json"
    CC_HOME="/tmp/cc_home_$$"
    mkdir -p "${CC_HOME}/.claude"

    # 创建最小 API 配置（DeepSeek via Anthropic 兼容接口），无会话历史
    # NOTE: credentials redacted for public release
    cat > "${CC_HOME}/.claude/settings.json" << 'CCJSON'
{
  "env": {
    "ANTHROPIC_AUTH_TOKEN": "<REDACTED:DEEPSEEK_API_KEY>",
    "ANTHROPIC_BASE_URL": "https://api.deepseek.com/anthropic",
    "ANTHROPIC_MODEL": "deepseek-v4-pro",
    "ANTHROPIC_SMALL_FAST_MODEL": "deepseek-v4-pro",
    "ANTHROPIC_DEFAULT_SONNET_MODEL": "deepseek-v4-pro",
    "ANTHROPIC_DEFAULT_OPUS_MODEL": "deepseek-v4-pro",
    "ANTHROPIC_DEFAULT_HAIKU_MODEL": "deepseek-v4-pro",
    "CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC": "1"
  }
}
CCJSON

    # bwrap 沙箱: 仅挂载 claudecode 工具 + 数据 + 工作区
    # 用环境变量传递 prompt（避免多行字符串在 bash -c 中断裂）
    BWRAP_ARGS=$(bwrap_common_args)
    bwrap \
      --ro-bind /workspace/claudecode /workspace/claudecode \
      ${BWRAP_ARGS} \
      --setenv HOME "${CC_HOME}" \
      --setenv PATH "/workspace/miniconda3/bin:/workspace/app/bin:/usr/bin:/bin:/workspace/claudecode/bin" \
      --setenv BENCH_PROMPT "${PROMPT}" \
      --setenv CC_AGENT_DIR "${AGENT_DIR}" \
      bash -c 'cd "$CC_AGENT_DIR" && /workspace/claudecode/bin/claude -p "$BENCH_PROMPT" --output-format json --allowedTools "Bash,Read,Write,Edit,Glob,Grep"' \
      > "${CC_RAW}" 2>&1
    EXIT_CODE=$?
    cp "${CC_RAW}" "${SESSION_LOG}"
    cp "${CC_RAW}" "${CONTEXT_JSON}"

    # Token from JSON
    if command -v jq >/dev/null 2>&1; then
      IN_TOK=$(jq -r '.usage.input_tokens // 0' "${CC_RAW}" 2>/dev/null || echo 0)
      OUT_TOK=$(jq -r '.usage.output_tokens // 0' "${CC_RAW}" 2>/dev/null || echo 0)
      CACHE_R=$(jq -r '.usage.cache_read_input_tokens // 0' "${CC_RAW}" 2>/dev/null || echo 0)
      CACHE_W=$(jq -r '.usage.cache_creation_input_tokens // 0' "${CC_RAW}" 2>/dev/null || echo 0)
      REAS_TOK=0
      SESSION_ID=$(jq -r '.session_id // ""' "${CC_RAW}" 2>/dev/null || echo "")
    fi

    # 清理临时 HOME
    mkdir -p "/tmp/ai_trash_$(date +%s)" && mv "${CC_HOME}" "/tmp/ai_trash_$(date +%s)/" 2>/dev/null || true
    ;;

  openclaw)
    # OpenClaw: bwrap 沙箱，无持久状态
    OC_RAW="${LOG_DIR}/oc_raw_$$.log"
    OC_HOME="/tmp/oc_home_$$"
    OC_SESSION="${TASK_ID}-$(date +%s)"
    mkdir -p "${OC_HOME}/.openclaw"

    # 复制完整配置到临时 HOME（含 model provider 定义，禁用 session-memory）
    # 同时拷贝 auth store（仅含 API key，非会话记忆）—— OpenClaw 凭此向 DeepSeek 鉴权
    mkdir -p "${OC_HOME}/.openclaw/agents/main/agent"
    if [ -f /root/.openclaw/agents/main/agent/openclaw-agent.sqlite ]; then
      cp /root/.openclaw/agents/main/agent/openclaw-agent.sqlite "${OC_HOME}/.openclaw/agents/main/agent/"
    fi
    # NOTE: gateway token redacted for public release
    cat > "${OC_HOME}/.openclaw/openclaw.json" << 'OCJSON'
{
  "agents": {
    "defaults": {
      "models": {
        "deepseek/deepseek-v4-pro": {"alias": "DeepSeek"}
      },
      "model": {"primary": "deepseek/deepseek-v4-pro"}
    }
  },
  "gateway": {
    "mode": "local",
    "auth": {"mode": "token", "token": "<REDACTED:LOCAL_GATEWAY_TOKEN>"},
    "port": 18789,
    "bind": "loopback"
  },
  "session": {"dmScope": "per-channel-peer"},
  "tools": {"profile": "coding"},
  "plugins": {"entries": {"deepseek": {"enabled": true}}},
  "auth": {"profiles": {"deepseek:default": {"provider": "deepseek", "mode": "api_key"}}},
  "models": {
    "mode": "merge",
    "providers": {
      "deepseek": {
        "baseUrl": "https://api.deepseek.com",
        "api": "openai-completions",
        "models": [
          {"id": "deepseek-v4-pro", "name": "DeepSeek V4 Pro", "reasoning": true, "input": ["text"], "contextWindow": 1000000, "maxTokens": 384000, "compat": {"supportsReasoningEffort": true, "supportsUsageInStreaming": true, "maxTokensField": "max_tokens"}, "api": "openai-completions"}
        ]
      }
    }
  },
  "hooks": {"internal": {"entries": {"session-memory": {"enabled": false}}}}
}
OCJSON

    # bwrap 沙箱: 仅挂载 openclaw 工具 + 数据 + 工作区
    BWRAP_ARGS=$(bwrap_common_args)
    # NOTE: DEEPSEEK_API_KEY must be provided via environment at run time
    bwrap \
      --ro-bind /workspace/openclaw /workspace/openclaw \
      ${BWRAP_ARGS} \
      --setenv HOME "${OC_HOME}" \
      --setenv PATH "/workspace/miniconda3/bin:/workspace/app/bin:/workspace/software/node-v24.18.0-linux-x64/bin:/usr/bin:/bin" \
      --setenv DEEPSEEK_API_KEY "${DEEPSEEK_API_KEY:-<REDACTED:DEEPSEEK_API_KEY>}" \
      --setenv BENCH_PROMPT "${PROMPT}" \
      --setenv OC_SESSION "${OC_SESSION}" \
      --setenv OC_AGENT_DIR "${AGENT_DIR}" \
      --setenv OC_PREFIX "/workspace/openclaw" \
      bash -c 'cd "$OC_AGENT_DIR" && npx --prefix "$OC_PREFIX" openclaw agent --agent main --session-key "$OC_SESSION" --local --json --message "$BENCH_PROMPT"' \
      > "${OC_RAW}" 2>&1
    EXIT_CODE=$?
    cp "${OC_RAW}" "${SESSION_LOG}"

    # 提取 JSON 响应
    OC_RAW_ENV="${OC_RAW}" python3 << 'PYEOF' > "${CONTEXT_JSON}" 2>/dev/null
import json, os, re
path = os.environ.get('OC_RAW_ENV')
with open(path) as f:
    text = f.read()
candidates = ['{"payloads"', '{"meta"', '{"sessionId"']
start = -1
for marker in candidates:
    idx = text.find(marker)
    if idx >= 0 and (start < 0 or idx < start):
        start = idx
if start < 0:
    for m in re.finditer(r'^\{', text, re.MULTILINE):
        if start < 0 or m.start() > start:
            start = m.start()
if start < 0:
    print('{"error": "no json start"}')
else:
    depth = 0
    in_string = False
    escape = False
    end = -1
    for i, ch in enumerate(text[start:], start):
        if escape:
            escape = False
            continue
        if ch == '\\':
            escape = True
            continue
        if ch == '"':
            in_string = not in_string
            continue
        if in_string:
            continue
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    if end > 0:
        try:
            data = json.loads(text[start:end])
            print(json.dumps(data))
        except json.JSONDecodeError as e:
            print(f'{{"error": "parse failed: {e}"}}')
    else:
        print('{"error": "no balanced json"}')
PYEOF
    [ -s "${CONTEXT_JSON}" ] || echo '{"error": "empty"}' > "${CONTEXT_JSON}"

    SESSION_ID="${OC_SESSION}"
    CTX_ENV="${CONTEXT_JSON}" python3 << 'PYEOF2' > "${LOG_DIR}/oc_token_$$.txt"
import json, os
with open(os.environ.get('CTX_ENV')) as f:
    try:
        data = json.load(f)
        usage = data.get('meta', {}).get('agentMeta', {}).get('usage', {})
        print(f"{usage.get('input', 0)} {usage.get('output', 0)} {usage.get('reasoningTokens', 0)} {usage.get('cacheRead', 0)}")
    except Exception:
        print('0 0 0 0')
PYEOF2
    TOKEN_INFO=$(cat "${LOG_DIR}/oc_token_$$.txt")
    read -r IN_TOK OUT_TOK REAS_TOK CACHE_R <<< "${TOKEN_INFO}"
    CACHE_W=0

    # 清理临时 HOME
    mkdir -p "/tmp/ai_trash_$(date +%s)" && mv "${OC_HOME}" "/tmp/ai_trash_$(date +%s)/" 2>/dev/null || true
    ;;

esac

END_TS=$(date +%s)
END_ISO=$(date -Iseconds)
DURATION=$((END_TS - START_TS))
TOTAL_TOK=$((IN_TOK + OUT_TOK + REAS_TOK))

# 收集产出文件清单
OUTPUT_FILES=$(find "${AGENT_DIR}" -type f -not -name "prompt.txt" -not -name "session_log.txt" -not -name "session_context.json" -not -name "metrics.json" 2>/dev/null | sed "s|${AGENT_DIR}/||" | sort | tr '\n' ',' | sed 's|,$||')

# 写 metrics.json
cat > "${METRICS_JSON}" <<EOF
{
  "task_id": "${TASK_ID}",
  "agent": "${AGENT}",
  "session_id": "${SESSION_ID}",
  "start_time": "${START_ISO}",
  "end_time": "${END_ISO}",
  "duration_sec": ${DURATION},
  "input_tokens": ${IN_TOK},
  "output_tokens": ${OUT_TOK},
  "cache_read_tokens": ${CACHE_R},
  "cache_write_tokens": ${CACHE_W},
  "reasoning_tokens": ${REAS_TOK},
  "total_tokens": ${TOTAL_TOK},
  "exit_code": ${EXIT_CODE},
  "completion_status": "pending",
  "output_files": [$(echo "${OUTPUT_FILES}" | tr ',' '\n' | sed 's|^|    "|; s|$|"|' | paste -sd, -)]
}
EOF

echo "[${TASK_ID}/${AGENT}] exit=${EXIT_CODE} duration=${DURATION}s in=${IN_TOK} out=${OUT_TOK} total=${TOTAL_TOK}"
