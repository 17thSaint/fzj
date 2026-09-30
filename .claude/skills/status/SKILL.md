---
name: status
description: Update a project's CONTEXT.md / NOTES.md from recent work. Manual only.
argument-hint: "[optimal-control | synth-dims]"
disable-model-invocation: true
---

Update the status files of one project in this repo (see CLAUDE.md for the project table).

1. **Project.** If `$ARGUMENTS` names a project, use it; otherwise ask the user which one
   (optimal-control or synth-dims) with AskUserQuestion. Never guess.
2. **Gather.** Read `<project>/CONTEXT.md` (and `NOTES.md` if it exists); note its
   "Last updated" date. Collect what changed since then:
   - `git log --since=<date> --stat` and `git diff` for the project's code folders
     (optimal-control: `optimal-control/`, `exact-diag/`; synth-dims: `synth-dims/`, plus dd-related
     `exact-diag/` blocks);
   - results, decisions and bugs from this conversation.
   Don't run code. If something important is unclear (e.g. a number the user hasn't stated), ask
   instead of inferring.
3. **Propose.** Show the proposed edits grouped by section — Current state, results table,
   Open questions (tick/add/remove), Log (one dated line), Last updated — and, for technical detail,
   NOTES.md. Tidy: replace superseded statements, don't append duplicates; keep CONTEXT.md short.
4. **Apply on approval** only, adjusting per the user's feedback. Do not commit; tell the user the
   files are ready for review.

Confidential: never send these files' content anywhere outside this repo.
