# CTOC14 Leaderboard Monitor

How the Problem A leaderboard evolved, rebuilt from the pages saved by hand under `leaderboard/<YYYYMMDD_HHMM>/`
(folder stamp = UTC; the board's times are Beijing). Each page lists the top 50 teams with rank, team, organisation,
time of the best submission, shown J and N.

## Rebuild the data

```
~/.venvs/astro313/bin/python tools/leaderboard_monitor.py --print
```

writes `board.json`, `board.js` (the same data as `window.BOARD`, read by the page so `file://` works) and
`events.csv` (one row per distinct submission). Re-run it after saving a new page.

## View

```
cd results/leaderboard && ~/.venvs/astro313/bin/python -m http.server 8766
```

then open <http://127.0.0.1:8766/>. Opening `index.html` directly from the file system also works (d3 is the only
external resource, from cdn.jsdelivr.net).

## What is shown

* **J chart**: shown J (what the board ranks by) or raw J (= shown / k, k = 1 - 0.1 (09-28 12:00 - t)/28 d) per team
  against submission time; lower is up; log or linear scale; each step is one submission (dot) held until the next.
  Line width = N (board N = spacecraft + missed targets; craft = N - 2 for every team at 298).
  "J range: current scores" fits the y axis to the current scores of the teams shown so the top-10 differences are visible
  (early history runs off the bottom); "full history" shows every step.
* **Rank chart**: rank by shown J, recomputed from the reconstructed histories at every submission time (exact between
  pages for every team a page ever listed).
* **Standings table**: current board; rows outside the "top" filter are dimmed; hover to highlight, click to pin.
* A team marked "gone" (dagger) disappeared from the page although its score would still rank (renamed, e.g. MSMLA ->
  MsMLab, or removed); it is ranked only until the last page that listed it and hidden unless "include removed teams".
* Left of the first saved page (shaded) only the entries later pages still showed are known.
