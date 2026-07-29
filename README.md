# PigGPA

Reproducible evidence and benchmark data for **PigGPA** (Genomic breeding; Prediction; Annotation), an AI agent system for pig genomics research. This repository contains three benchmark subsets:

- **agent_benchmark/** — PigGPA vs. general-purpose coding agents (Claude Code, OpenClaw) on a real GWAS + regulatory-activity-prediction task (G→P→A workflow), scored on a 10-point rubric.
- **piggpa_benchmark/** — End-to-end smoke test of all 26 PigGPA-G sub-skills (100% pass), 100-query intent parsing (100% pass), and 5 error-handling scenarios (100% pass).
- **piggpa-hiblup_benchmark/** — Quantitative comparison of piggpa-G vs. HIBLUP (v1.6.0) across 15 genomic-prediction modules, with metric-by-metric verification on identical simulated pig population data.
