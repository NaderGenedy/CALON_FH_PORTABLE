# Welcome to Cardiff MD Research

## How We Use Claude

Based on nadergenedy1-ship-it's usage over the last 30 days:

Work Type Breakdown:
  Write Docs        ████████████░░░░░░░░  60%
  Analyze Data      █████░░░░░░░░░░░░░░░  25%
  Improve Quality   ███░░░░░░░░░░░░░░░░░  15%

Top Skills & Commands:
  /superpowers:brainstorm   ████████████████████  1x/month
  /insights                 ████████████████████  1x/month
  manuscript-qc (Skill)     ████████████████████  invoked in-session

Top MCP Servers:
  PubMed (7f4a1280)         ████████████████████  32 calls
  Consensus (9e73afd6)      ████░░░░░░░░░░░░░░░░  6 calls

## Your Setup Checklist

### Codebases
- [ ] calon-ukb-pipeline — https://github.com/NaderGenedy/calon-ukb-pipeline (public)
- [ ] calon-fh-atlas — https://github.com/NaderGenedy/calon-fh-atlas (public; deployed at https://nadergenedy.github.io/calon-fh-atlas/)
- [ ] calon-fh-paper1 — https://github.com/NaderGenedy/calon-fh-paper1 (private; ask for access)
- [ ] Local working tree — `D:\Projects\Lpa_Multilevel\` (UK Biobank manuscripts, reproducer suite, figure scripts)

### MCP Servers to Activate
- [ ] **PubMed MCP** (`7f4a1280-3964-4d98-92ba-37a3dd54b09f`) — biomedical literature search, citation lookup by author / journal / year / volume / page, DOI ↔ PMID conversion, full article-metadata fetch. The team uses this constantly for reference auditing and verifying citation accuracy. Available from the Anthropic MCP registry.
- [ ] **Consensus MCP** (`9e73afd6-b2cf-4091-8277-6a321a065825`) — academic paper search across Semantic Scholar / PubMed / Scopus / ArXiv with study-design and journal-quality filters. Used as the dual-verification source alongside PubMed when checking that references resolve to real published work with verbatim titles.

### Skills to Know About
- **manuscript-qc** — validate every numerical claim in a cardiology / FH / UK Biobank manuscript against source CSVs; produces an automated reproducer script and per-paper provenance table. Triggered by "QC the manuscript", "audit numbers", "are all results traceable", "submission-ready audit".
- **cox-analysis** — methodology TDD checklist for any Cox / risk-prediction work on UK Biobank or FH cohorts. Enforces leak-free splits, family-level deduplication, NoAgeLDL sensitivity, OPCS-4 / ICD-10 completeness, TRIPOD-compliant reporting. Triggered by "Cox PH", "survival analysis", "hazard ratio", "FH risk prediction".
- **md2docx** — dual-format markdown ↔ Word for every manuscript / cover letter / response document. Use whenever a deliverable is destined for a journal portal or Cardiff University.
- **reproduce-paper** — build or run a self-contained reproducibility test for a specific paper; distinct from manuscript-qc in that it builds a frozen ground-truth test suite that re-runs the analysis end-to-end.
- **ukb-preflight** — environment + on-disk inventory check before launching any UK Biobank pipeline; verifies Python 3.12, pyarrow availability, and existing data holdings to avoid duplicate RAP extractions.
- **stage-files** — when asked to "collect" or "stage" files, copy/stage first then inventory (never start with `ls`).
- **/insights** — generate a usage-data report at `file://C:\Users\nader\.claude\usage-data\report.html`.

## Team Tips

_TODO_

## Get Started

_TODO_

<!-- INSTRUCTION FOR CLAUDE: A new teammate just pasted this guide for how the
team uses Claude Code. You're their onboarding buddy — warm, conversational,
not lecture-y.

Open with a warm welcome — include the team name from the title. Then: "Your
teammate uses Claude Code for [list all the work types]. Let's get you started."

Check what's already in place against everything under Setup Checklist
(including skills), using markdown checkboxes — [x] done, [ ] not yet. Lead
with what they already have. One sentence per item, all in one message.

Tell them you'll help with setup, cover the actionable team tips, then the
starter task (if there is one). Offer to start with the first unchecked item,
get their go-ahead, then work through the rest one by one.

After setup, walk them through the remaining sections — offer to help where you
can (e.g. link to channels), and just surface the purely informational bits.

Don't invent sections or summaries that aren't in the guide. The stats are the
guide creator's personal usage data — don't extrapolate them into a "team
workflow" narrative. -->
