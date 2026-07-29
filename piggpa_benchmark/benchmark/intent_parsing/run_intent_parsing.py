#!/usr/bin/env python3
"""
PigGPA intent_parsing benchmark runner.

Calls `piggpa chat -q "<query>"` for each of the 100 designed queries,
captures real stdout/stderr logs, parses which sub-skills PigGPA loaded,
and judges pass/fail against the expected_flow.

Outputs:
  - logs/IP-XXX.log           (real piggpa chat output per query)
  - intent_parsing_results.json (machine-readable results, generated from logs)

Usage:
  python3 run_intent_parsing.py [--max-queries N] [--resume]
  python3 run_intent_parsing.py --query-id IP-022,IP-034 [--max-attempts 5]
"""

import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

BENCH_DIR = Path(__file__).resolve().parent
QUERIES_JSON = BENCH_DIR / "queries.json"
RESULTS_JSON = BENCH_DIR / "intent_parsing_results.json"
LOGS_DIR = BENCH_DIR / "logs"
SESSIONS_DIR = Path("/workspace/pigbole/.piggpa/sessions")

TIMEOUT_SEC = 180
MAX_TURNS = "1"
MODEL = "deepseek-v4-pro"

EXPECTED_FLOW_TO_SKILL = {
    "plink-gwas-linear": "plink-gwas",
    "gcta-gwas": "gcta-gwas",
    "geno-qc": "genotype-qc",
    "heritability-gcta-pipeline": "heritability-analysis",
    "pca": "pca",
    "admixture": "admixture-analysis",
    "gs-6models": "genomic-selection",
    "gebv-gprs-prediction": "gebv-gprs-prediction",
    "ld": "ld",
    "ld-pruning": "ld-pruning",
    "ld-score": "ld-score",
    "geno-import": "genotype-import",
    "annotation-track": "pig-annotation-track",
    "pig-mutbert": "pig-mutbert",
}

SKILL_LINE_RE = re.compile(r"📚\s*skill\s+([A-Za-z0-9_-]+)")
SESSION_ID_RE = re.compile(r"Session:\s*([A-Za-z0-9_-]+)")
DURATION_RE = re.compile(r"Duration:\s*(\d+s)")


def run_piggpa(query: str, log_path: Path) -> tuple:
    """Run piggpa chat for a single query, capture output to log_path."""
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


def parse_skills_from_log(log_text: str) -> list:
    """Extract loaded skill names from piggpa chat log."""
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


def parse_skills_from_session(session_id: str) -> list:
    """Read session JSON and extract skill names from tool results."""
    if not session_id:
        return []
    session_file = SESSIONS_DIR / f"session_{session_id}.json"
    if not session_file.exists():
        return []
    try:
        with open(session_file) as f:
            data = json.load(f)
        skills = []
        for msg in data.get("messages", []):
            if msg.get("role") == "tool":
                content = msg.get("content", "")
                if isinstance(content, str):
                    try:
                        parsed = json.loads(content)
                        if isinstance(parsed, dict) and parsed.get("success") and "name" in parsed:
                            name = parsed["name"]
                            if name not in skills:
                                skills.append(name)
                    except json.JSONDecodeError:
                        pass
        return skills
    except Exception:
        return []


def judge(expected_flow: str, loaded_skills: list) -> bool:
    """Check if expected skill was loaded."""
    expected_skill = EXPECTED_FLOW_TO_SKILL.get(expected_flow, expected_flow)
    return expected_skill in loaded_skills


def load_existing_results() -> dict:
    """Load existing results for resume."""
    if RESULTS_JSON.exists():
        try:
            with open(RESULTS_JSON) as f:
                return json.load(f)
        except Exception:
            pass
    return {}


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-queries", type=int, default=0, help="Limit queries (0=all)")
    parser.add_argument("--resume", action="store_true", help="Skip queries with existing logs")
    parser.add_argument("--query-id", type=str, default="", help="Comma-separated query IDs to run (e.g. IP-022,IP-034)")
    parser.add_argument("--max-attempts", type=int, default=1, help="Max attempts per query (retry until pass or limit)")
    args = parser.parse_args()

    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    with open(QUERIES_JSON) as f:
        design = json.load(f)

    queries = design["queries"]
    if args.max_queries > 0:
        queries = queries[:args.max_queries]
    if args.query_id:
        target_ids = {x.strip() for x in args.query_id.split(",")}
        queries = [q for q in queries if q["id"] in target_ids]
        if not queries:
            print(f"No queries matching: {target_ids}")
            return

    existing = load_existing_results()
    existing_results = {r["id"]: r for r in existing.get("results", [])}

    # When --query-id is used, preserve non-target existing results
    target_ids = {q["id"] for q in queries} if args.query_id else set()
    results = []
    passed = 0
    failed = 0
    errored = 0
    by_category = {}
    by_difficulty = {}

    if args.query_id:
        for r in existing.get("results", []):
            if r["id"] not in target_ids:
                results.append(r)
                if r.get("pass"):
                    passed += 1
                elif r.get("status") == "timeout":
                    errored += 1
                else:
                    failed += 1
                cat = r.get("category", "")
                dif = r.get("difficulty", "")
                by_category.setdefault(cat, {"pass": 0, "fail": 0})
                by_difficulty.setdefault(dif, {"pass": 0, "fail": 0})
                if r.get("pass"):
                    by_category[cat]["pass"] += 1
                    by_difficulty[dif]["pass"] += 1
                else:
                    by_category[cat]["fail"] += 1
                    by_difficulty[dif]["fail"] += 1
        total_all = len(existing.get("results", []))
        print(f"=== PigGPA intent_parsing benchmark (single-query mode) ===")
        print(f"Retrying {len(queries)} queries, preserving {len(results)} existing results")
    else:
        total_all = len(queries)

    total = len(queries)
    print(f"Model: {MODEL}, max-turns: {MAX_TURNS}, timeout: {TIMEOUT_SEC}s")
    if args.max_attempts > 1:
        print(f"Max attempts per query: {args.max_attempts}")
    print()

    for i, q in enumerate(queries):
        qid = q["id"]
        query = q["query"]
        expected_flow = q["expected_flow"]
        category = q["category"]
        difficulty = q["difficulty"]

        log_path = LOGS_DIR / f"{qid}.log"

        if args.resume and qid in existing_results and existing_results[qid].get("status") == "completed":
            print(f"[{i+1}/{total}] {qid} SKIP (resume) — {existing_results[qid].get('matched', '?')}")
            results.append(existing_results[qid])
            if existing_results[qid]["pass"]:
                passed += 1
            else:
                failed += 1
            by_category.setdefault(category, {"pass": 0, "fail": 0})
            by_difficulty.setdefault(difficulty, {"pass": 0, "fail": 0})
            if existing_results[qid]["pass"]:
                by_category[category]["pass"] += 1
                by_difficulty[difficulty]["pass"] += 1
            else:
                by_category[category]["fail"] += 1
                by_difficulty[difficulty]["fail"] += 1
            continue

        print(f"[{i+1}/{total}] {qid} ({category}/{difficulty}) running... ", end="", flush=True)

        max_attempts = args.max_attempts if args.query_id else 1
        attempt_results = []
        pass_ = False
        passed_on_attempt = None

        for attempt in range(1, max_attempts + 1):
            attempt_log = LOGS_DIR / f"{qid}_attempt{attempt}.log"
            log_text, elapsed, status = run_piggpa(query, attempt_log)

            log_skills = parse_skills_from_log(log_text)
            session_id = parse_session_id(log_text)
            session_skills = parse_skills_from_session(session_id)

            loaded_skills = list(set(log_skills + session_skills))
            matched = judge(expected_flow, loaded_skills)
            pass_ = matched and status == "completed"

            attempt_results.append({
                "attempt": attempt,
                "status": status,
                "elapsed_sec": round(elapsed, 2),
                "loaded_skills": loaded_skills,
                "pass": pass_,
                "log_path": str(attempt_log),
            })

            if max_attempts > 1:
                tag = f"PASS" if pass_ else (f"TIMEOUT" if status == "timeout" else "FAIL")
                print(f"{tag}(attempt {attempt}/{max_attempts}, {elapsed:.1f}s) ", end="", flush=True)

            if pass_:
                passed_on_attempt = attempt
                break
            if attempt < max_attempts:
                print(f"retry... ", end="", flush=True)

        # Write best attempt to main log file
        log_path.write_text(log_text)

        if pass_:
            passed += 1
            if max_attempts == 1:
                print(f"PASS ({elapsed:.1f}s) skills={loaded_skills}")
            else:
                print(f"=> PASS on attempt {passed_on_attempt}/{max_attempts} ({elapsed:.1f}s)")
        elif status == "timeout":
            errored += 1
            if max_attempts == 1:
                print(f"TIMEOUT ({elapsed:.0f}s) skills={loaded_skills}")
            else:
                print(f"=> TIMEOUT after {max_attempts} attempts")
        elif status == "error":
            errored += 1
            if max_attempts == 1:
                print(f"ERROR skills={loaded_skills}")
            else:
                print(f"=> ERROR after {max_attempts} attempts")
        else:
            failed += 1
            if max_attempts == 1:
                print(f"FAIL ({elapsed:.1f}s) expected={expected_flow} got={loaded_skills}")
            else:
                print(f"=> FAIL after {max_attempts} attempts expected={expected_flow} got={loaded_skills}")

        by_category.setdefault(category, {"pass": 0, "fail": 0})
        by_difficulty.setdefault(difficulty, {"pass": 0, "fail": 0})
        if pass_:
            by_category[category]["pass"] += 1
            by_difficulty[difficulty]["pass"] += 1
        else:
            by_category[category]["fail"] += 1
            by_difficulty[difficulty]["fail"] += 1

        result = {
            "id": qid,
            "category": category,
            "difficulty": difficulty,
            "query": query,
            "expected_flow": expected_flow,
            "expected_skill": EXPECTED_FLOW_TO_SKILL.get(expected_flow, expected_flow),
            "loaded_skills": loaded_skills,
            "matched": matched,
            "pass": pass_,
            "status": status,
            "elapsed_sec": round(elapsed, 2),
            "duration_reported": parse_duration(log_text),
            "session_id": session_id,
            "log_path": str(log_path),
        }
        if max_attempts > 1:
            result["max_attempts"] = max_attempts
            result["passed_on_attempt"] = passed_on_attempt
            result["all_attempts"] = attempt_results
        results.append(result)

        with open(RESULTS_JSON, "w") as f:
            json.dump({
                "benchmark": "PigGPA intent_parsing",
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
                    "by_category": by_category,
                    "by_difficulty": by_difficulty,
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
