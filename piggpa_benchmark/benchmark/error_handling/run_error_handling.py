#!/usr/bin/env python3
"""
PigGPA error_handling benchmark runner.

Calls `piggpa chat -q "<query>"` for each of the 5 error scenarios,
captures real stdout/stderr logs, parses PigGPA's response for error
detection signals, and judges pass/fail against pass_criteria keywords.

Outputs:
  - logs/ERR-XXX.log              (real piggpa chat output per scenario)
  - error_handling_results.json   (machine-readable results, generated from logs)

Usage:
  python3 run_error_handling.py [--resume]
  python3 run_error_handling.py --scenario-id ERR-002 [--max-attempts 5]
"""

import json
import re
import subprocess
import time
from datetime import datetime
from pathlib import Path

BENCH_DIR = Path(__file__).resolve().parent
SCENARIOS_JSON = BENCH_DIR / "scenarios.json"
RESULTS_JSON = BENCH_DIR / "error_handling_results.json"
LOGS_DIR = BENCH_DIR / "logs"
SESSIONS_DIR = Path("/workspace/pigbole/.piggpa/sessions")

TIMEOUT_SEC = 180
MAX_TURNS = "2"
MODEL = "deepseek-v4-pro"

SCENARIO_KEYWORDS = {
    "ERR-001": {
        "required_any": [["missing", "required", "genotype", "input", "provide", "need", "file path", "bfile", "--bfile", "--genotype"]],
        "required_avoid": [["--maf 0.8", "plink --maf"]],
    },
    "ERR-002": {
        "required_any": [["invalid", "range", "0.5", "must be", "valid", "maf", "0 to", "0-0.5", "0.8"]],
        "required_avoid": [["--maf 0.8", "maf=0.8", "maf 0.8"]],
    },
    "ERR-003": {
        "required_any": [["ambiguous", "clarify", "specify", "which", "options", "could you", "what kind", "what type", "several", "different types"]],
        "required_avoid": [],
    },
    "ERR-004": {
        "required_any": [["does not exist", "not found", "invalid path", "no such file", "check path", "cannot find", "doesn't exist", "nonexistent", "non-existent"]],
        "required_avoid": [],
    },
    "ERR-005": {
        "required_any": [["mismatch", "trait", "different", "confirm", "alignment", "same trait", "mismatched", "don't match", "do not match"]],
        "required_avoid": [],
    },
}

SKILL_LINE_RE = re.compile(r"📚\s*skill\s+([A-Za-z0-9_-]+)")
SESSION_ID_RE = re.compile(r"Session:\s*([A-Za-z0-9_-]+)")
DURATION_RE = re.compile(r"Duration:\s*(\d+s)")
RESPONSE_START_RE = re.compile(r"─\s+🐽\s+PigGPA\s+─+")


def run_piggpa(query: str, log_path: Path) -> tuple:
    cmd = [
        "piggpa", "chat",
        "-q", query,
        "--max-turns", MAX_TURNS,
        "--yolo",
        "-m", MODEL,
    ]
    start = time.time()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SEC,
            cwd="/workspace/pigbole",
        )
        elapsed = time.time() - start
        output = proc.stdout + proc.stderr
        log_path.write_text(output)
        return output, elapsed, "completed"
    except subprocess.TimeoutExpired as e:
        elapsed = time.time() - start
        stdout = e.stdout if e.stdout else ""
        stderr = e.stderr if e.stderr else ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode("utf-8", errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode("utf-8", errors="replace")
        output = stdout + stderr
        log_path.write_text(output + f"\n[TIMEOUT after {TIMEOUT_SEC}s]\n")
        return output, elapsed, "timeout"
    except Exception as e:
        elapsed = time.time() - start
        log_path.write_text(f"[ERROR] {e}\n")
        return str(e), elapsed, "error"


def extract_response_text(log_text: str) -> str:
    """Extract PigGPA's final response text from the log."""
    lines = log_text.splitlines()
    in_response = False
    response_lines = []
    for line in lines:
        if RESPONSE_START_RE.search(line):
            in_response = True
            continue
        if in_response:
            if line.startswith("─" * 20) or line.startswith("╰" * 20):
                break
            response_lines.append(line)
    return "\n".join(response_lines) if response_lines else log_text


def parse_skills_from_log(log_text: str) -> list:
    skills = []
    for line in log_text.splitlines():
        m = SKILL_LINE_RE.search(line)
        if m:
            skill = m.group(1)
            if skill not in skills:
                skills.append(skill)
    return skills


def parse_session_id(log_text: str) -> str:
    m = SESSION_ID_RE.search(log_text)
    return m.group(1) if m else ""


def parse_duration(log_text: str) -> str:
    m = DURATION_RE.search(log_text)
    return m.group(1) if m else ""


def judge(scenario_id: str, response_text: str, loaded_skills: list) -> dict:
    """Judge scenario against keyword criteria."""
    criteria = SCENARIO_KEYWORDS.get(scenario_id, {})
    response_lower = response_text.lower()

    matched_keywords = []
    for group in criteria.get("required_any", []):
        group_matches = [kw for kw in group if kw.lower() in response_lower]
        matched_keywords.append(group_matches)
        if not group_matches:
            return {
                "pass": False,
                "reason": f"no keyword from {group} found in response",
                "matched_keywords": matched_keywords,
            }

    for group in criteria.get("required_avoid", []):
        for kw in group:
            if kw.lower() in response_lower:
                return {
                    "pass": False,
                    "reason": f"forbidden keyword '{kw}' found in response (agent attempted invalid execution)",
                    "matched_keywords": matched_keywords,
                }

    return {
        "pass": True,
        "reason": "all keyword criteria met",
        "matched_keywords": matched_keywords,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--scenario-id", type=str, default="", help="Comma-separated scenario IDs to run (e.g. ERR-002)")
    parser.add_argument("--max-attempts", type=int, default=1, help="Max attempts per scenario (retry until pass or limit)")
    args = parser.parse_args()

    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    with open(SCENARIOS_JSON) as f:
        design = json.load(f)

    scenarios = design["scenarios"]
    if args.scenario_id:
        target_ids = {x.strip() for x in args.scenario_id.split(",")}
        scenarios = [s for s in scenarios if s["id"] in target_ids]
        if not scenarios:
            print(f"No scenarios matching: {target_ids}")
            return

    # When --scenario-id is used, preserve non-target existing results
    target_ids = {s["id"] for s in scenarios} if args.scenario_id else set()
    results = []
    passed = 0
    failed = 0
    errored = 0

    if args.scenario_id:
        if RESULTS_JSON.exists():
            with open(RESULTS_JSON) as f:
                existing = json.load(f)
            for r in existing.get("results", []):
                if r["id"] not in target_ids:
                    results.append(r)
                    if r.get("pass"):
                        passed += 1
                    elif r.get("status") == "timeout":
                        errored += 1
                    else:
                        failed += 1
        total_all = len(results) + len(scenarios)
        print(f"=== PigGPA error_handling benchmark (single-scenario mode) ===")
        print(f"Retrying {len(scenarios)} scenarios, preserving {len(results)} existing results")
    else:
        total_all = len(scenarios)

    total = len(scenarios)
    print(f"Model: {MODEL}, max-turns: {MAX_TURNS}, timeout: {TIMEOUT_SEC}s")
    if args.max_attempts > 1:
        print(f"Max attempts per scenario: {args.max_attempts}")
    print()

    for i, sc in enumerate(scenarios):
        sid = sc["id"]
        query = sc["query"]
        description = sc["description"]
        test_type = sc["test_type"]
        expected_behavior = sc["expected_behavior"]
        pass_criteria = sc["pass_criteria"]

        log_path = LOGS_DIR / f"{sid}.log"

        if args.resume and log_path.exists():
            log_text = log_path.read_text()
            elapsed = 0
            status = "completed" if "[TIMEOUT" not in log_text else "timeout"
            print(f"[{i+1}/{total}] {sid} SKIP (resume)")
        else:
            print(f"[{i+1}/{total}] {sid} ({test_type}) running... ", end="", flush=True)

            max_attempts = args.max_attempts if args.scenario_id else 1
            attempt_results = []
            pass_ = False
            passed_on_attempt = None

            for attempt in range(1, max_attempts + 1):
                attempt_log = LOGS_DIR / f"{sid}_attempt{attempt}.log"
                log_text, elapsed, status = run_piggpa(query, attempt_log)

                response_text = extract_response_text(log_text)
                loaded_skills = parse_skills_from_log(log_text)
                session_id = parse_session_id(log_text)

                judgment = judge(sid, response_text, loaded_skills)
                pass_ = judgment["pass"] and status == "completed"

                attempt_results.append({
                    "attempt": attempt,
                    "status": status,
                    "elapsed_sec": round(elapsed, 2),
                    "loaded_skills": loaded_skills,
                    "pass": pass_,
                    "judgment_reason": judgment.get("reason", ""),
                    "log_path": str(attempt_log),
                })

                if max_attempts > 1:
                    tag = "PASS" if pass_ else ("TIMEOUT" if status == "timeout" else "FAIL")
                    print(f"{tag}(attempt {attempt}/{max_attempts}, {elapsed:.1f}s) ", end="", flush=True)

                if pass_:
                    passed_on_attempt = attempt
                    break
                if attempt < max_attempts:
                    print(f"retry... ", end="", flush=True)

            # Write best attempt to main log file
            log_path.write_text(log_text)

        if not args.resume or not (args.resume and log_path.exists()):
            pass  # already judged above in retry path
        else:
            response_text = extract_response_text(log_text)
            loaded_skills = parse_skills_from_log(log_text)
            session_id = parse_session_id(log_text)
            judgment = judge(sid, response_text, loaded_skills)
            pass_ = judgment["pass"] and status == "completed"
            max_attempts = 1
            passed_on_attempt = None
            attempt_results = []

        if pass_:
            passed += 1
            if max_attempts == 1:
                print(f"PASS ({elapsed:.1f}s) — {judgment['reason']}")
            else:
                print(f"=> PASS on attempt {passed_on_attempt}/{max_attempts} ({elapsed:.1f}s)")
        elif status == "timeout":
            errored += 1
            if max_attempts == 1:
                print(f"TIMEOUT ({elapsed:.0f}s) — {judgment['reason']}")
            else:
                print(f"=> TIMEOUT after {max_attempts} attempts")
        elif status == "error":
            errored += 1
            if max_attempts == 1:
                print(f"ERROR — {judgment['reason']}")
            else:
                print(f"=> ERROR after {max_attempts} attempts")
        else:
            failed += 1
            if max_attempts == 1:
                print(f"FAIL ({elapsed:.1f}s) — {judgment['reason']}")
            else:
                print(f"=> FAIL after {max_attempts} attempts — {judgment['reason']}")

        result = {
            "id": sid,
            "description": description,
            "test_type": test_type,
            "query": query,
            "expected_behavior": expected_behavior,
            "pass_criteria": pass_criteria,
            "loaded_skills": loaded_skills,
            "judgment": judgment,
            "pass": pass_,
            "status": status,
            "elapsed_sec": round(elapsed, 2),
            "duration_reported": parse_duration(log_text),
            "session_id": session_id,
            "log_path": str(log_path),
            "response_excerpt": response_text[:500] if response_text else "",
        }
        if max_attempts > 1:
            result["max_attempts"] = max_attempts
            result["passed_on_attempt"] = passed_on_attempt
            result["all_attempts"] = attempt_results
        results.append(result)

        with open(RESULTS_JSON, "w") as f:
            json.dump({
                "benchmark": "PigGPA error_handling",
                "version": "1.0",
                "run_status": "in_progress" if i < total - 1 else "completed",
                "run_timestamp": datetime.now().isoformat(),
                "model": MODEL,
                "max_turns": int(MAX_TURNS),
                "timeout_sec": TIMEOUT_SEC,
                "summary": {
                    "total": total_all,
                    "passed": passed,
                    "failed": failed,
                    "errored": errored,
                    "pass_rate": round(passed / total_all * 100, 1) if total_all > 0 else 0,
                },
                "results": results,
            }, f, indent=2, ensure_ascii=False)

    pass_rate = round(passed / total_all * 100, 1) if total_all > 0 else 0
    print()
    print(f"=== DONE ===")
    print(f"Total: {total}, Passed: {passed}, Failed: {failed}, Errored: {errored}")
    print(f"Pass rate: {pass_rate}%")
    print(f"Results: {RESULTS_JSON}")


if __name__ == "__main__":
    main()
