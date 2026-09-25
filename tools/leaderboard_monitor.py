"""Leaderboard monitor: parse the saved educoder pages under leaderboard/<YYYYMMDD_HHMM>/ (UTC stamp in the folder
name, board times are Beijing = UTC+8) and write results/leaderboard/{board.json, board.js, events.csv}.
The page results/leaderboard/index.html reads board.js (inline data, so file:// works; only d3 comes from a CDN).

    ~/.venvs/astro313/bin/python tools/leaderboard_monitor.py            # rebuild the data
    ~/.venvs/astro313/bin/python tools/leaderboard_monitor.py --print    # also print the current standings

Every board row is (rank, team, org, time of the best submission, shown J, N).  The shown J of a submission never
changes afterwards (checked over all snapshots), so a team's history is the set of distinct (submit time, J, N) rows,
and its shown score is a step function of time.  raw J = shown / k(t_submit), k = 1 - 0.1 (t_end - t)/28 d.
Board N = craft + 2 for every team that reaches 298 targets (2 are unreachable): it is the number of spacecraft plus
the number of missed targets.  Ranks are recomputed from the reconstructed histories at every submission time
(the pages only rank at snapshot times), so the rank track is exact between snapshots for every team the pages ever
listed (top 50 only, so a team can be missing before it first entered the top 50)."""
import re, csv, html, json, glob, pathlib, argparse, datetime as dt

ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
SRC = ROOT / 'leaderboard'
OUT = ROOT / 'results/leaderboard'
T_START = dt.datetime(2026, 8, 31, 12, 0)
T_END = dt.datetime(2026, 9, 28, 12, 0)
US = 'Mickey'


def k_of(t):
    return 1.0 - 0.1 * (T_END - t).total_seconds() / (T_END - T_START).total_seconds()


def parse_page(path):
    s = open(path, encoding='utf-8', errors='replace').read()
    rows = []
    for tr in re.findall(r'<tr class="ant-table-row[^"]*">(.*?)</tr>', s, re.S):
        tds = re.findall(r'<td[^>]*>(.*?)</td>', tr, re.S)
        if len(tds) < 6:
            continue
        name = re.sub(r'<div[^>]*>.*?</div>', '', tds[1], flags=re.S)       # drop the avatar initial
        cells = [html.unescape(re.sub(r'<[^>]+>', '', x)).strip() for x in [tds[0], name] + tds[2:6]]
        rows.append(dict(rank=int(cells[0]), team=cells[1], org=cells[2], sub=cells[3], J=float(cells[4]), N=int(cells[5])))
    return rows


def iso(t):
    return t.strftime('%Y-%m-%d %H:%M')


def build():
    snaps = []
    for d in sorted(SRC.glob('2026*_*')):
        f = d / 'ctoc14.educoder.net.html'
        if not f.exists():
            continue
        t = dt.datetime.strptime(d.name, '%Y%m%d_%H%M') + dt.timedelta(hours=8)     # folder stamp is UTC
        snaps.append(dict(id=d.name, t=iso(t), rows=parse_page(f)))
    teams = {}
    for si, sn in enumerate(snaps):
        for r in sn['rows']:
            tm = teams.setdefault(r['team'], dict(team=r['team'], org=r['org'], events={}, first_snap=si, last_snap=si))
            tm['last_snap'] = si
            key = (r['sub'], r['J'], r['N'])
            ev = tm['events'].setdefault(key, dict(t=r['sub'], J=r['J'], N=r['N'], first_seen=sn['t'], seen=0))
            ev['seen'] += 1
            ev['last_seen'] = sn['t']
            ev['rank_seen'] = r['rank']
    # a submission time seen with two different J values would break the step model: report it
    for tm in teams.values():
        by_t = {}
        for (sub, J, N) in tm['events']:
            by_t.setdefault(sub, set()).add((J, N))
        for sub, v in by_t.items():
            if len(v) > 1:
                print(f'WARNING {tm["team"]}: submission {sub} shown with several values {sorted(v)}')
    out_teams = []
    for tm in teams.values():
        evs = sorted(tm['events'].values(), key=lambda e: e['t'])
        for e in evs:
            t = dt.datetime.strptime(e['t'], '%Y-%m-%d %H:%M')
            e['k'] = round(k_of(t), 6)
            e['raw'] = round(e['J'] / e['k'], 6)
        last = evs[-1]
        out_teams.append(dict(team=tm['team'], org=tm['org'], events=evs, n_sub=len(evs),
                              first_snap=snaps[tm['first_snap']]['t'], last_snap=snaps[tm['last_snap']]['t'],
                              dropped=tm['last_snap'] != len(snaps) - 1,
                              cur=dict(t=last['t'], J=last['J'], raw=last['raw'], N=last['N'], k=last['k'])))
    # rank track: at every submission time, order all teams that have submitted by then (by shown J, the board's key).
    # A team that vanished from the page although its J would still rank (renamed, e.g. MSMLA -> MsMLab, or removed)
    # is ranked only until the last snapshot that listed it.
    times = sorted({e['t'] for tm in out_teams for e in tm['events']} | {s['t'] for s in snaps})
    tracks = []
    for t in times:
        cur = []
        for tm in out_teams:
            if tm['dropped'] and t > tm['last_snap']:
                continue
            evs = [e for e in tm['events'] if e['t'] <= t]
            if evs:
                cur.append((evs[-1]['J'], tm['team']))
        cur.sort()
        tracks.append(dict(t=t, order=[name for _, name in cur]))
    out_teams.sort(key=lambda tm: (tm['dropped'], tm['cur']['J']))
    for i, tm in enumerate(out_teams, 1):
        tm['cur']['rank'] = i if not tm['dropped'] else None
    data = dict(generated=iso(dt.datetime.now()), t_start=iso(T_START), t_end=iso(T_END), us=US,
                n_snapshots=len(snaps), first_snapshot=snaps[0]['t'], last_snapshot=snaps[-1]['t'],
                snapshot_times=[s['t'] for s in snaps], teams=out_teams, rank_track=tracks,
                snapshots=[dict(id=s['id'], t=s['t'], rows=s['rows']) for s in snaps])
    return data


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--print', action='store_true'); a = ap.parse_args()
    data = build()
    OUT.mkdir(parents=True, exist_ok=True)
    json.dump(data, open(OUT / 'board.json', 'w'), ensure_ascii=False, indent=0)
    with open(OUT / 'board.js', 'w') as f:
        f.write('window.BOARD = ' + json.dumps(data, ensure_ascii=False) + ';\n')
    with open(OUT / 'events.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['team', 'org', 'submit_time', 'J_shown', 'k', 'J_raw', 'N_board', 'first_seen', 'last_seen', 'rank_when_seen'])
        for tm in data['teams']:
            for e in tm['events']:
                w.writerow([tm['team'], tm['org'], e['t'], e['J'], e['k'], e['raw'], e['N'], e['first_seen'], e['last_seen'], e['rank_seen']])
    n_ev = sum(tm['n_sub'] for tm in data['teams'])
    print(f'{data["n_snapshots"]} snapshots {data["first_snapshot"]} .. {data["last_snapshot"]}, {len(data["teams"])} teams, '
          f'{n_ev} distinct submissions -> {OUT}/board.js')
    if a.print:
        print(f'{"rk":>3} {"team":<22} {"J shown":>10} {"J raw":>9} {"N":>3}  {"last submit":<16} subs')
        for tm in data['teams'][:20]:
            c = tm['cur']
            print(f'{c["rank"]:>3} {tm["team"]:<22} {c["J"]:>10.6f} {c["raw"]:>9.4f} {c["N"]:>3}  {c["t"]:<16} {tm["n_sub"]}'
                  + ('  (dropped from page)' if tm['dropped'] else ''))


if __name__ == '__main__':
    main()
