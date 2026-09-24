# CTOC / CTOC14 contest background

Research snapshot taken 2026-09-06 (UTC+8 dates unless stated). Sources are listed inline and in §9.
Everything about CTOC14 itself was taken from (i) the official CSTAM notice of 2026-08-03 and (ii) the JSON
API behind the competition site `ctoc14.educoder.net`; the problem PDF in this repo remains authoritative for the
problem itself. §8 lists what could NOT be found.

---

## 0. Executive summary for the lead

| Item | Value | Source |
|---|---|---|
| Competition | 第十四届中国空间轨道设计竞赛 (CTOC14), 14th China Trajectory Optimization Competition | [N1], [E1] |
| Host | 主办: 中国力学学会 (CSTAM). 承办: 国防科技大学空天科学学院, 太空系统运行与控制全国重点实验室 (NUDT), 空天技术国家级实验教学示范中心, 清华大学行健书院, 清华大学航天航空学院 | [N1], [E3-28] |
| Problem A (甲题) setter | 太空系统运行与控制全国重点实验室 (State Key Lab of Space System Operation and Control, NUDT, Changsha). Lab personnel may not enter 甲题. 甲题负责人: 黄岸毅 (email hay04@qq.com) | [N1], [E3-30], [E3-33] |
| Problem A theme | 近地小行星防御多航天器遍历探测轨道设计与优化 (multi-spacecraft low-thrust flyby tour of 300 NEAs; notice wording: "面向近地小行星探测任务的轨道设计与优化") | PDF, [N1] |
| Problem B (乙题) | 故障卫星巡检轨道设计与优化 (inspection of a failed satellite), undergraduate-only, solved in NUDT's ATK software, hosted by Tsinghua; opens 2026-09-21 12:00 | [N1], [E2] |
| **t_start (甲题)** | **2026-08-31 12:00 Beijing time (UTC+8) = 2026-08-31 04:00 UTC** | [E2], [E3-29], [N1] |
| **t_end (甲题)** | **2026-09-28 12:00 Beijing time = 2026-09-28 04:00 UTC** (also "竞赛结束, 公布结果") | [E2], [E3-29], [N1] |
| Duration | exactly 28.0 days = 672 h | derived |
| Registration | opens 2026-08-01 12:00, closes 2026-08-30 12:00 per notice; the site's own record says enrol 2026-07-24 17:31 → 2026-09-28 00:00 (site data overrides notice in practice, but do not rely on late registration) | [N1], [E2] |
| Website | https://ctoc14.educoder.net (EduCoder/头歌 platform instance; SPA, content served from `/api/competitions/...`) | [N1], [E1] |
| Submission | 甲题: one `.txt` file (`CTOC14_Result_TeamID.txt`), uploaded on the website; **online validation limited to 30 runs per team per day (12:00 → next-day 12:00)**; e-mail submission only as a fallback if the site crashes | [E3-31] |
| Time coefficient | k = 1 − 0.1·(t_end − t_submit)/(t_end − t_start); k = 0.9 at t_start, 1.0 at t_end; **computed from the LAST valid submission** (an early submission does not "lock in" a low k if you resubmit later) | PDF §5.3 |
| Prizes 甲题 | 1st 特等奖 ¥15,000; 2–3 一等奖 ¥6,000; 4–6 二等奖 ¥4,000; 7–10 三等奖 ¥1,500; 11–20 优胜奖 certificate. 甲题 champion earns the right to host 甲题 of CTOC15 | [E3-32], [N1] |
| Team rules 甲题 | no limit on team size or member level; each supervisor may supervise at most one 甲题 team; one account per team | [E3-30] |
| Validation | organizer's independent checker re-integrates every line with 3rd-order sliding-window Lagrange-interpolated thrust; tolerances 1 km / 1 m/s / 0.01 kg; ‖T‖ ≤ 0.5 N checked continuously; the online checker on the site is the same code family used for 乙题/ATK (for 乙题 the statement is explicit; for 甲题 the PDF describes "独立程序") | PDF §7, [E3-31] |
| Site activity | 158 registered member accounts and ~14,000 page visits as of 2026-09-06 | [E1] |

**Practical consequence of k (details in §5):** each day of delay raises the final score by 0.1/28 = 0.357 % of J
(absolute: J/280 per day). Keep optimizing only while you expect to cut J by more than ≈0.36–0.39 % per day
(e.g. for J ≈ 100 that is ≈ 0.36–0.39 cost units per day, roughly "one extra asteroid every ~2.6 days" or
"one fewer spacecraft every ~2.6–7.8 days" depending on its fuel load).

---

## 1. What CTOC is

* **Name.** 中国空间轨道设计竞赛 / 全国空间轨道设计竞赛, "China Trajectory Optimization Competition" (CTOC).
  Founded 2009 as 全国深空轨道设计竞赛 ("national deep-space trajectory design competition"); renamed 全国空间轨道设计竞赛
  from the 4th edition (2012) when problems stopped being limited to deep space; since ~2023 official material also
  uses 中国空间轨道设计竞赛. [R10]
* **Organizer.** 中国力学学会 (Chinese Society of Theoretical and Applied Mechanics, CSTAM) is the permanent 主办单位;
  the competition "依托中国力学学会". Official portal: https://ctoc.cstam.org.cn/ (archive pages CTOC-1…CTOC-12,
  统章程/组委会; it is a JS-rendered site — see §8). [R10], [C0]
* **Origin.** Proposed in 2009 by Prof. 李俊峰 (Li Junfeng) and Prof. 宝音贺西 (Baoyin Hexi) of Tsinghua's School of Aerospace
  Engineering, explicitly modelled on ESA's GTOC (2005–). CSTAM's 刘俊丽 picked it up in April 2009 and CSTAM + Tsinghua
  co-hosted the first edition. [R10]
* **Rotation rule.** As in GTOC, **the champion team's institution sets and hosts the next edition** (for two-track editions:
  the 甲题 champion hosts the next 甲题; since CTOC13 the undergraduate-track champion hosts the next undergraduate track). [R10], [N13], [N1]
* **Cadence.** Annual 2009–2017; every 1–2 years afterwards (2019, 2020, 2022, 2024, 2026). [R10]
* **Format.** Open problem with an objective, machine-checkable performance index; no reference answer (organizers admit
  they cannot find the global optimum either). Solving window 30–60 days (CTOC9: 45 d; CTOC11: 50 d; CTOC12: 31 d for A;
  CTOC13: 31 d; CTOC14: 28 d for 甲题), versus GTOC's 28 days. Teams may use any human/computational resources. Organizer
  re-computes every submitted trajectory and ranks purely by the index. Since CTOC13 there is a live online leaderboard and
  online self-validation. [R10], [R12], [U13]
* **Standing.** The 10-year review credits CTOC with raising Chinese teams from the bottom half of GTOC1–4 to winning GTOC10
  (2019: NUDT + 西安卫星测控中心 first, Tsinghua + 星邑空间 second). Results are traditionally written up in 力学与实践
  (Mechanics in Engineering), with CTOC9 in a special section of *Acta Astronautica* and CTOC10 in *Astrodynamics*. [R10]

## 2. Edition-by-edition history (CTOC1 – CTOC13)

Sources: 10-year review [R10] (CTOC1–10), CTOC11 news [U11], CTOC12 review [R12]/[B12]/[S12], CTOC13 notice [N13],
CTOC13-B winners' paper [P13B], CTOC13 undergraduate summary [U13], symposium news [SY13], CTOC14 notice [N1].

| # | Year (problem release → deadline) | Host / setter (besides CSTAM) | Problem(s) | Teams (registered / valid) | Champion(s) |
|---|---|---|---|---|---|
| 1 | 2009 (2009-03-05 → 04-16) | 清华大学 (Tsinghua) | 小行星探测轨迹优化设计 (asteroid exploration / sample-return trajectory) | 33 / – | 中国科学院光电研究院 (CAS Academy of Opto-Electronics); 2nd 西安卫星测控中心, 3rd 国防科技大学 |
| 2 | 2010 (03-25 → ~May) | 西安卫星测控中心 + 中科院空间科学与应用总体部 (then part of 光电研究院) | Mars + multiple NEA multi-target exploration | 28 / 12 (15 submitted) | 清华大学航天航空学院 (winner used electric propulsion and planetary gravity assists) |
| 3 | 2011 (03-10 →) | 北京航天飞行控制中心 + 航天飞行动力学技术重点实验室 + Tsinghua | Tour of the eight planets with small-body visits incl. asteroid "钱学森" (Qian Xuesen centenary) | 25 / 11 | 中国科学院空间应用工程与技术中心 (CAS CSU) |
| 4 | 2012 (06-10 →) | CAS CSU + BACC lab + NUDT 航天与材料工程学院 | Multi-target (asteroids + comets) multi-task (flyby, rendezvous, impact, sample return) small-body exploration; renamed 全国空间轨道设计竞赛 | 26 / 13 (14 submitted) | 国防科技大学 (NUDT; 沈红新, 罗亚中, 李海阳) |
| 5 | 2013 (06-01 →) | NUDT 航天科学与工程学院 + 航天飞行动力学技术重点实验室 | Manned near-Earth-asteroid exploration trajectory | 20 / 13 | tie: 西安卫星测控中心宇航动力学国家重点实验室, 中科院光电研究院, CAS CSU |
| 6 | 2014 (07-15 →) | 西安卫星测控中心宇航动力学国家重点实验室 | First two-track edition. A: NEA sample return in multi-body gravity field; B: fastest escape from the solar system | – | A: CAS CSU; B: tie NUDT (朱阅訸, 罗亚中, 贺波勇) and 光电研究院+CSU alliance |
| 7 | 2015 (08-15 →) | CAS CSU | A: rover tour on the surface of an irregular asteroid; B: LEO satellite formation reconfiguration | 40 / 19 | "一等奖" (no single champion): A tie NSSC, 西安 lab, Tsinghua; B tie NSSC, 光电研究院, 西安 lab |
| 8 | 2016 (08-08 →) | 中科院国家空间科学中心 (NSSC) + 西安 lab | A: multi-target SSO space-debris removal; B: satellite planning/scheduling for multi-target ground observation | 62 / 22 | A: 北京理工大学/航天二院二部 alliance; B: tie Tsinghua 宇航中心 and CAS CSU |
| 9 | 2017 (09-01 → 10-15, 45 d) | BIT/航天二院二部 + Tsinghua 宇航中心 | A: GEO satellite beam-monitoring cluster flight; B: regional navigation-augmentation constellation design/deployment. First bilingual edition (3 foreign teams incl. ESA) | 52 / 18 (21 submitted) | A: 南京航空航天大学 (NUAA); B: NSSC. Papers in Acta Astronautica 150 (2018) |
| 10 | 2019 (03-20 → 04-30) | NUAA 航天学院 | Jupiter internal magnetic field + Galilean moon science exploration (multi-probe, low thrust, gravity assists) | 62 (8 foreign, 4 commercial) / 20 | 哈尔滨工业大学 (HIT). First commercial sponsor (北京宇航智科). Summary: Astrodynamics 5(1):1–11 (2021) |
| 11 | 2020 (50 days, results 2020-11-09) | HIT | Dual-satellite coordinated orbit design for large-scale ground-target tasks: 0–60 d full coverage of 200 static targets; 60–180 d coverage + equal-interval revisit; 180–240 d repeated uniform observation of 5 moving targets | 49 (incl. Jena Univ.) / 9 | "空间应用中心队" (CAS CSU alumni + UCAS students, 何胜茂 lead), score 602 |
| 12 | 2022 (reg. from 05-20; A: 08-20 → 09-20; B: 08-20 → 10-10) | 中科院太空应用重点实验室 (CAS CSU) + 国科大航空宇航学院 | A: reusable Earth–Mars transport stations based at lunar DRO carrying migrants to Mars (≤50 stations, 20 yr); **B: reusable probes from cislunar space flying by hazardous NEAs, max asteroids in 10 yr** | 39 / 10 per problem | A: 北京理工大学-航天东方红卫星有限公司-北京工业大学 alliance (9,080 migrants, 50 stations, 519 km/s); 2nd Σ团队 (NSSC + BACC) 5,540/28; 3rd 北航–航天科工二院 4,320/24. **B: Σ团队 (NSSC + BACC): 45 asteroids with 11 probes, 1.80 km/s per asteroid**; 2nd 南京大学–国家天文台 39/11; 3rd 36/16; 4th 四川大学 15/6; 5th NUAA 9/4 |
| 13 | 2024 (reg. 07-01 12:00; A/B released 08-20; C released 09-13; deadline 09-20 12:00) | BIT + 航天东方红卫星有限公司 + 北京工业大学 + CAS CSU + NUDT (A set by BIT 温昶煊, B by CSU 何胜茂, C by NUDT 朱阅訸) | A: 巨型星座巡检轨道设计与优化 (multi-flyby inspection of a mega-constellation); B: 地月可重复使用运输飞行器轨道设计 (LEO ↔ lunar-DRO reusable tankers, MJD 60676–61406, ≤20 vehicles); C (first undergraduate track, solved in ATK): 面向应急对地观测任务的混合星座设计与机动调度 | C alone: 79 teams from 19 universities | A: **not confirmed** (see §8; inferred NUDT / 太空系统运行与控制全国重点实验室 because the A champion hosts CTOC14 甲题); B: 北京航天飞行控制中心 (BACC; 李皓皓, 刘勇 et al.); C: Tsinghua "力3biu" (two 2nd-year undergraduates), 2nd/3rd NUDT. Prizes A ¥20k/10k/5k, B ¥10k/5k/3k. Symposium + 1st ATK user conference, Changsha, 2025-04-12/13, 270+ attendees |

Notes on the table:
* Through CTOC10 the review counts "14 problems, champions from 10 institutions". CAS CSU's 空间探索室 claims 6 championships
  (2009, 2011, 2013, 2014-A, 2016-B, 2020) and 3 hostings (2010, 2012, 2015). [U11]
* Institutions that recur as winners/hosts and are therefore the likely strongest CTOC14 competitors: Tsinghua (Li Junfeng /
  Baoyin Hexi / Jiang Fanghua group), NUDT (Luo Yazhong / Yang Zhen / Zhu Yuehe — excluded from 甲题 this year as hosts),
  CAS CSU (Gao Yang / He Shengmao), NSSC (asteroid-defence group), BACC (Liu Yong), 西安卫星测控中心, BIT (Zhang Jingrui),
  NUAA (Li Shuang / Yang Hongwei), HIT, BUAA, 航天东方红.
* CTOC12-B (reusable NEA flyby probes, 45 asteroids in 10 yr with 11 probes) and CTOC1/2/4/5 are the closest CTOC ancestors
  of CTOC14 甲题. CTOC12-B winners relied on lunar gravity assists for inclination changes, constrained-v∞ transfers,
  consecutive-flyby optimization and "distributed resonant flyby touring orbits"; CTOC14 forbids gravity assists, so only the
  latter ideas carry over. [R12]

## 3. CTOC14 (2026) — everything public

### 3.1 Announcement
* CSTAM notice "第十四届中国空间轨道设计竞赛通知", published 2026-08-03 at https://m.cstam.org.cn/186/202608/22039.html [N1].
  The competition site record was created 2024-09-05, went online 2026-06-30 16:21 and was last edited 2026-08-27 [E1].
* Stated purpose: "为促进学术交流，进一步推动空间轨道设计技术的发展，更好服务于我国未来的空间探索事业，计划于 2026年8月31日至9月28日举办第十四届中国空间轨道设计竞赛。"

### 3.2 Organizers
* 主办: 中国力学学会.
* 承办: 国防科技大学空天科学学院; 太空系统运行与控制全国重点实验室 (NUDT, Changsha 410073 — the lab that also co-hosted the CTOC13
  symposium and whose staff 朱阅訸/杨震/罗亚中 ran CTOC13's undergraduate track and develop ATK); 空天技术国家级实验教学示范中心;
  清华大学行健书院 (dean 李俊峰); 清华大学航天航空学院. [N1], [SY13], [U13]
* 甲题承办单位 = 太空系统运行与控制全国重点实验室 (its personnel barred from 甲题); 乙题承办单位 = 清华大学 (its students barred from 乙题). [E3-30]
* The EduCoder tenant record lists `sponsor_schools: ["火箭军工程大学"]` — an artefact of the hosting platform account, not
  an organizer named in the notice. [E1]

### 3.3 Themes (notice wording)
* 甲题: "小行星撞击是地球面临的现实潜在威胁之一……本届竞赛甲题即设定为：面向近地小行星探测任务的轨道设计与优化。" The PDF title is
  近地小行星防御多航天器遍历探测轨道设计与优化. [N1], PDF
* 乙题: "剧烈太阳风暴和空间碎片撞击等对在轨高价值卫星的安全运行构成了严重威胁……本届竞赛乙题设定为：故障卫星巡检轨道设计与优化。"
  Undergraduates only (incl. 2026 graduates), ≤5 per team, solved with a competition build of ATK released with the problem. [N1], [E3-30]

### 3.4 Schedule (日程安排, Beijing time UTC+8)
| Event | Date/time |
|---|---|
| 注册报名开启 | 2026-08-01 12:00 |
| 注册报名截止 | 2026-08-30 12:00 |
| **甲题竞赛开始，发布赛题 (t_start)** | **2026-08-31 12:00** |
| 乙题竞赛开始，发布赛题 | 2026-09-21 12:00 |
| **竞赛结束，公布结果 (t_end)** | **2026-09-28 12:00** |

The site's topic record for 甲题 carries identical `start_time: 2026-08-31 12:00`, `end_time: 2026-09-28 12:00` [E2]; the
PDF says t_start/t_end are "以竞赛网站公布为准", so these are the values to use. Results are announced at t_end.

### 3.5 Registration and team rules (注册报名)
1. 甲题: no restriction on number or level of members. 乙题: undergraduates only, ≤5 members.
2. One member registers the team account; in "账号管理中心" choose 甲题/乙题 and fill team name, institution, supervisor, members.
   **Each supervisor may supervise at most one 甲题 team and one 乙题 team.** A single account may enter both tracks. [E3-30]

### 3.6 Submission and validation (结果提交)
1. 甲题 result = `.txt` file strictly in the PDF's 15-column format; 乙题 result = `.atk` scene file. For 乙题 the site's
   back-end validator "与ATK软件中的验证程序为同一套代码".
2. **Submissions must go through the website; each team may run online validation at most 30 times per 24 h window
   (12:00 → next 12:00).**
3. If the site crashes, results may be e-mailed and will be validated/scored manually. [E3-31]
4. PDF §7: checker verifies format/structure, dynamics consistency (re-integration line-to-line, ≤1 km, ≤1 m/s, ≤0.01 kg),
   launch (position ≤1 km from Earth, v∞ ≤ 4 km/s with 0.01 km/s tolerance), thrust (‖T‖ ≤ 0.5 N at every instant of the
   interpolated profile), detection (d ≤ 1000 km at declared Event=3 lines; first detection only), mass (m0 ≤ 2000, m ≥ 600),
   and computes N_covered, N_miss = 300 − N_covered and J.
5. Ranking: ascending J_final = k·J; "任务总代价 J 越小、提交越早，排名越靠前". Multiple submissions allowed; k uses the last valid one;
   submissions after t_end are ignored. PDF §5.3–5.4

### 3.7 Prizes (评奖办法)
Each track ranks and awards independently, 20 awarded teams per track:

| Rank | Level | 甲题 | 乙题 |
|---|---|---|---|
| 1 | 特等奖 | ¥15,000 | ¥10,000 |
| 2–3 | 一等奖 | ¥6,000 | ¥5,000 |
| 4–6 | 二等奖 | ¥4,000 | ¥3,000 |
| 7–10 | 三等奖 | ¥1,500 | ¥1,000 |
| 11–20 | 优胜奖 | certificate | certificate |

The champion of each track gains the right to host that track of CTOC15. [E3-32], [N1]

### 3.8 Contacts and reference material
* 甲题负责人 黄岸毅 (hay04@qq.com); 乙题负责人 张楠; competition WeChat group via 朱阅訸 (NUDT). Phone numbers are in the notice [N1]/[E3-33].
* Reference downloads offered for newcomers (undergraduate track only): CTOC13 丙题 package and the 15th 周培源力学竞赛
  space-trajectory team problem, both with ATK builds (OSS links in [N1]). Latest ATK: https://www.osredm.com/atknudt/atk/about.
* Site modules "赛题背景" and "竞赛题目" were empty at fetch time (problem delivered as the PDF + MEA.txt after 2026-08-31 12:00); the
  "竞赛通知" module contains only "……". No informs/announcements have been posted on the site (informs.json = []). Leaderboard
  (`charts.json`, `chart_rules.json`) requires login. [E3-34/35/36]

## 4. CTOC13 → CTOC14 lineage and what it implies

* CTOC13's notice stated that the 甲题 champion would host CTOC14 [N13]; CTOC14 甲题 is hosted by NUDT's 太空系统运行与控制全国重点实验室,
  so that lab (or a team from it) is very probably the CTOC13 甲题 (mega-constellation inspection) champion — unconfirmed (§8).
  The lab's public CTOC record: NUDT won CTOC4 (2012, multi-target small-body tour, 沈红新/罗亚中/李海阳), CTOC6-B (2014, 朱阅訸/罗亚中),
  and was 2nd/3rd in CTOC13-C. Its research staff (罗亚中 Luo Yazhong, 杨震 Yang Zhen, 朱阅訸 Zhu Yuehe, 沈红新 Shen Hongxin) publish on
  evolutionary/global optimization of multi-target missions, rendezvous planning, and the ATK mission-analysis toolkit. [R10], [U13], [SY13]
* Undergraduate track lineage: CTOC13-C champion Tsinghua "力3biu" → Tsinghua hosts CTOC14 乙题. [U13], [N1]
* The 10-year review notes setters typically spend ~6 months designing and "反复试算" the problem to remove loopholes, and that
  problems are designed so that they cannot be solved with off-the-shelf tools. Expect the 300-asteroid catalogue to have been
  pre-screened so that the baseline (a few spacecraft, tens of asteroids each) is feasible but full coverage is hard. [R10]

## 5. Time coefficient k — numbers for decision making

k = 1 − 0.1·(t_end − t_submit)/(t_end − t_start), with t_end − t_start = 28 d exactly; k rises by 0.1/28 = 0.0035714 per day
(0.000149 per hour). J_final = k·J, evaluated at the **last valid** submission.

| Submit at (Beijing 12:00) | days elapsed | k |
|---|---|---|
| 2026-09-01 | 1 | 0.9036 |
| 2026-09-06 | 6 | 0.9214 |
| 2026-09-10 | 10 | 0.9357 |
| 2026-09-14 | 14 | 0.9500 |
| 2026-09-17 | 17 | 0.9607 |
| 2026-09-21 | 21 | 0.9750 |
| 2026-09-24 | 24 | 0.9857 |
| 2026-09-27 | 27 | 0.9964 |
| 2026-09-28 12:00 | 28 | 1.0000 |

* Delaying the final submission by Δt days is worth it only if J drops by more than a factor k(t)/k(t+Δt), i.e. by more than
  100·(Δt/280)/(k+Δt/280) % ≈ **0.39 % per day at k≈0.92, 0.36 % per day at k≈1.0**. In cost units: ≈ J/280 per day
  (J = 80 → 0.29/day; J = 120 → 0.43/day). One asteroid = 1.0 unit; an empty spacecraft = 1.0 unit; a full-tank spacecraft = 3.0 units.
* Because k is taken from the last valid submission, there is no benefit in submitting early except as insurance against an
  invalid/late final upload (an invalid later upload leaves the earlier valid one, with its earlier k, as the scored entry — the PDF says
  "以最后一次有效提交的时刻计算"). Submitting an improved file always resets k to the new (later) time.
* The 30-validations-per-day cap on the site means the local validator must be trusted; do not burn online validations on
  incremental changes.

## 6. Typical tooling and methods in CTOC (for problem-A-type multi-target low-thrust tours)

Sources: winners' write-ups and summaries [R10], [R12], [P13B], [B12], [U11], [S12], the GTOC/CTOC methodology review [PAS18],
and the 2025 Tsinghua/PoliMi global-optimality paper that uses GTOC4 as a benchmark [ZZ25]. (The full text of [PAS18] could not be
retrieved — see §8 — so its content is summarized from the abstract/title and from the papers that cite it.)

1. **Two-level decomposition** (universal in CTOC and GTOC): an outer combinatorial search over target sequences/timing on a cheap
   surrogate, and an inner continuous optimization of each low-thrust leg. Winners of CTOC13-B describe exactly this: "建立了LEO至DRO往返
   轨道数据库，通过整数优化求解多飞行器重复往返最优飞行序列" (transfer database + integer programming for multi-vehicle sequencing), with
   differential evolution (DE) for the 13-variable leg optimization and a "时间与轨道协同同伦" (time–orbit homotopy) continuation to widen
   convergence. [P13B]
2. **Surrogates for legs:** Lambert/bi-impulse Δv and pork-chop grids; orbital-element-based Δv estimates (Edelbaum-type,
   |Δa|,|Δe|,|Δi| metrics); for flyby-only legs a position-match (no velocity constraint) Lambert with limited Δv or a "reachability"
   test. The GTOC4 winner (MSU) first solved an **impulsive simplification with manoeuvres only at flyby points** (48 flybys + 1 rendezvous
   in that model) and then converted to low thrust (44 + 1 finally). Tsinghua's 2025 work replaces surrogates with neural-network
   estimators of low-thrust Δv/reachability and a reduced-dimension dynamic-programming search with provable error bounds. [ZZ25], [GT4]
3. **Global sequence search:** beam search / tree search with pruning by phase and Δv, dynamic programming over (time, target)
   grids, ant-colony and genetic algorithms, branch-and-bound, and integer/assignment models for multi-vehicle allocation
   (CTOC12-A/B and CTOC13-B used explicit "database + optimization" approaches; CTOC11's winner used a bespoke technical route with
   a very high first-phase score then global refinement over 50 days). [R12], [P13B], [U11]
4. **Local low-thrust optimization:** indirect methods (PMP, bang-bang fuel-optimal control with homotopy from energy-optimal —
   Tsinghua's signature tool), direct transcription (Sims–Flanagan segments or collocation solved with SQP/IPOPT), and hybrid shooting +
   differential correction (BUAA in CTOC12-A: "天体借力、弱稳定边界等手段，结合数值搜索、优化与微分修正"). [B12], [ZZ25]
5. **Multi-vehicle allocation:** greedy assignment then exchange/local search; treating the number of vehicles as an outer variable;
   CTOC12-B's winner used 11 probes for 45 asteroids and Nanjing Univ. 11 for 39 — i.e. 3–4 asteroids per probe in a
   cislunar-return setting with impulsive Δv; CTOC14's no-return low-thrust spacecraft should do better per spacecraft. [R12]
6. **Software:** overwhelmingly self-written MATLAB / C++ / Fortran propagators and optimizers; NUDT's ATK (Aerospace Tool Kit,
   v3.5 released April 2025, with a built-in result verifier identical to the competition back-end for the ATK tracks); the official
   checker's rules (Lagrange 3rd-order sliding-window thrust interpolation, 8640 s minimum sample spacing) must be replicated locally
   because the site allows only 30 online checks/day. [U13], [SY13], [E3-31]
7. **Historical yardsticks for "asteroids per spacecraft":** GTOC4 (one probe, 1435-NEA catalogue, 0.135 N, 1500 kg, 10 yr) → 44
   flybys + rendezvous in competition, 49 + 1 best known; CTOC12-B (impulsive, cislunar-based, 10 yr) → 45 asteroids with 11 probes.
   CTOC14's 300-target catalogue is ~5× sparser than GTOC4's, but its spacecraft have ~2.8× the acceleration (0.25 vs 0.09 mm/s² at
   full mass), higher Isp (4000 vs 3000 s) and 15 instead of 10 years, and a 1000 km detection sphere instead of exact position match.

## 7. GTOC4 (2009) — the closest international analogue, briefly

* **Organizer:** CNES (Bertrand Regis / Richard Epenoy), released 2009-03-02. Portal: https://sophia.estec.esa.int/gtoc_portal/?page_id=23 [GT4].
* **Problem:** one low-thrust spacecraft launched from Earth (v∞ ≤ 4 km/s), fly by as many distinct NEAs as possible (flyby = position match,
  any relative velocity) and end with a rendezvous at one NEA; ≤ 10 years from launch to rendezvous; launch window 2015–2025.
  Parameters (Table 5 of [ZZ25], quoting the problem statement): μ = 1.32712440018e20 m³/s², Isp = 3000 s, g0 = 9.80665, Tmax = 0.135 N,
  m0 = 1500 kg, minimum final mass 500 kg, ΔV_E ≤ 4 km/s. Catalogue: 1435 asteroids (`gtoc4_problem_data.txt`, epoch MJD 54800, J2000
  ecliptic Keplerian elements) [GT4D]. Objective: number of asteroids; tie-breakers final mass, then flight time.
* **Results:** 1st Moscow State University (Ilia S. Grigoriev) 44 flybys + 1 rendezvous; 2nd The Aerospace Corporation, also 44 but lower
  final mass. Post-competition: UT Austin 45+1; Jena University 46+1 (2015), 47+1 (2017), 49+1 (2018); [ZZ25] re-optimized fuel on the 49+1
  sequence (19.9 kg propellant left). [GT4], [ALT], [ZZ25]
* **Method of the winner:** impulsive-at-flyby simplification for sequence search, then low-thrust conversion; ranking PDFs on the portal
  are image/Type-3 encoded and could not be text-extracted (see §8). [ZZ25]
* Another agent covers GTOC4 methods in depth; the key transferable lesson is the two-level (sequence surrogate → low-thrust refinement)
  workflow and the fact that flyby-only tours sustain ~4–5 targets per year for a single low-thrust probe in a dense catalogue.

## 8. What could NOT be found / verified

* **CTOC13 甲题 champion** (mega-constellation inspection): not confirmed. `ctoc13.spacestdc.com` (the CTOC13 site with the leaderboards)
  no longer resolves; the symposium report [SY13] does not name winners; only the inference from the hosting rule is available.
* **Full text of the official CTOC portal** (https://ctoc.cstam.org.cn/, `/history`, `/competition-regulations`, `/news/230`, `/news/240–251`):
  the site is a JS-rendered "rhhz" journal-platform SPA; curl returns only navigation, the WebFetch tool refuses the domain, and the JSON
  endpoint could not be discovered. The 竞赛章程 (formal regulations) therefore remain unread; the CTOC13 notice text was taken from the
  thepaper.cn mirror [N13].
* **CTOC14 problem A team count / current leaderboard:** `charts.json`, `competition_teams.json`, `results.json` require login. Only
  `member_count = 158` (accounts) and `visits ≈ 14,000` are public.
* **CTOC14 乙题 details beyond the theme** (it opens 2026-09-21).
* **Whether ATK's validator is literally the same code as the 甲题 online checker** — stated only for 乙题 (.atk files); the PDF says the
  organizer uses "独立程序" for 甲题.
* **[PAS18] full text** (ScienceDirect 403; Crossref has no abstract) — only bibliographic data confirmed: Li Shuang, Huang Xuxing, Yang Bin,
  "Review of optimization methodologies in global and China trajectory optimization competitions", Progress in Aerospace Sciences 102
  (2018) 60–75, DOI 10.1016/j.paerosci.2018.07.004.
* **GTOC4 official ranking table and problem-description PDFs**: downloaded but not text-extractable (fonts without ToUnicode); values above
  come from the portal page, [ZZ25] and the data file.
* **CTOC1–CTOC10 exact team counts for some editions and the CTOC11 problem-setter paper** (a 2025 力学与实践 summary of CTOC11 was hinted
  at by search snippets but no URL/DOI was recovered).
* Web search quota was exhausted mid-task (session limit), so later items were obtained only via direct fetches of known URLs.

## 9. Sources

Official CTOC14
* [N1] 中国力学学会, 第十四届中国空间轨道设计竞赛通知, 2026-08-03. https://m.cstam.org.cn/186/202608/22039.html
* [E1] CTOC14 site record (JSON): https://ctoc14.educoder.net/api/competitions.json
* [E2] CTOC14 header/topics with per-track times (JSON): https://ctoc14.educoder.net/api/competitions/h6xcbmrw/common_header.json
* [E3-nn] CTOC14 site modules (JSON, nn = 28 竞赛简介, 29 日程安排, 30 注册报名, 31 结果提交, 32 评奖办法, 33 联系我们, 34 赛题背景, 35 竞赛题目, 36 竞赛通知):
  https://ctoc14.educoder.net/api/competitions/h6xcbmrw/competition_modules/28.json … /36.json (site: https://ctoc14.educoder.net/)
* PDF: /Users/mickey/solarsystem/ctoc14/CTOC14_problem.pdf (第十四届全国空间轨道设计竞赛（CTOC14）题目甲, 太空系统运行与控制全国重点实验室)

CTOC history
* [C0] 中国空间轨道设计竞赛 official portal (CSTAM): https://ctoc.cstam.org.cn/ ; 竞赛介绍 https://ctoc.cstam.org.cn/competition-introduction ;
  CTOC13 notice https://ctoc.cstam.org.cn/news/230 ; CTOC11 page https://ctoc.cstam.org.cn/news/250
* [R10] 刘俊丽, 高扬. 十年一剑刃锋利, 苦寒方得梅花香—全国空间轨道设计竞赛发展历程回顾. 力学与实践, 2019, 41(4): 488–497. DOI 10.6052/1000-0879-19-283.
  https://pubs.cstam.org.cn/article/doi/10.6052/1000-0879-19-283?viewType=HTML (full HTML also at
  http://lxsj.cstam.org.cn/fileLXYSJ/journal/article/lxysj/html/article/2019/1000-0879/1000-0879-41-4-488.shtml). Its reference list gives the
  per-edition summaries in 力学与实践: CTOC2 2011,33(2):116; CTOC3 2012,34(2):97 and 34(3):95; CTOC4 2013,35(1):99 and 2012,34(6):95;
  CTOC5 2014,36(3):379; CTOC6 2015,37(2):276/282 and 37(4):557; CTOC7 2016,38(5):596 and 38(6):99; CTOC9 Acta Astronautica 150 (2018) 177–249.
* [3C] 清华大学/国防科技大学/中科院空间应用中心, 创建"3C模式"专业竞赛，培养航天创新人才 成果总结报告 (2016):
  https://www.hy.tsinghua.edu.cn/__local/5/FF/F9/0508EAE24386164B0E28FA6E25F_57481324_7666F.pdf
* [A10] Huang X., Yang B., Sun P., Li S., Yang H. 10th China Trajectory Optimization Competition: Problem description and summary of the results.
  Astrodynamics 5(1):1–11 (2021). DOI 10.1007/s42064-020-0089-2. https://www.sciopen.com/article/10.1007/s42064-020-0089-2
* [U11] 中国科学院大学星际航行学院, 航空宇航学院研究生荣获第十一届全国空间轨道设计竞赛冠军 (2020-11-17):
  https://sse.ucas.ac.cn/index.php/zh/xyxw/120-2020-11-17-04-25-17
* [R12] 第十二届全国（中国）空间轨道设计竞赛：地月空间助力火星移民和近地小行星探测. 力学与实践 (2024), DOI 10.6052/1000-0879-23-603.
  https://pubs.cstam.org.cn/article/doi/10.6052/1000-0879-23-603?viewType=HTML (mirror: https://lihang.dlut.edu.cn/info/1070/16151.htm)
* [B12] 北航宇航学院, 我院师生荣获全国空间轨道设计竞赛季军 (2022-09-22): http://www.sa.buaa.edu.cn/info/1050/8822.htm
* [S12] 四川大学空天科学与工程学院, CTOC12 乙组第四名: https://saa.scu.edu.cn/cdkt/content.jsp?id=661626750534551
* [N13] 关于举办第十三届空间轨道设计竞赛的通知 (mirror of the CSTAM notice): https://www.thepaper.cn/newsDetail_forward_27901254
* [P13B] 李皓皓, 刘勇, 李革非, 等. 第十三届中国空间轨道设计竞赛乙题冠军团队解法. 力学与实践, 2026, 48(3): 585–594. DOI 10.6052/1000-0879-25-243.
  https://lxsj.cstam.org.cn/en/article/pdf/preview/10.6052/1000-0879-25-243.pdf
* [U13] 朱阅訸, 杨震, 罗亚中. 中国空间轨道设计竞赛首次本科生赛道举办探索情况总结. 力学与实践, 2025, 47(3): 684–692. DOI 10.6052/1000-0879-25-039.
  https://pubs.cstam.org.cn/article/doi/10.6052/1000-0879-25-039
* [SY13] 中国力学学会, 第十三届中国空间轨道设计竞赛研讨会暨首届ATK软件用户大会顺利召开 (2025-04-25): https://m.cstam.org.cn/187/202504/21334.html
* [CA13] arXiv 2507.02943, Global Optimization of Multi-Flyby Trajectories for Multi-Orbital-Plane Constellations Inspection (probably a CTOC13-A write-up; not read)

Methods / GTOC
* [PAS18] Li S., Huang X., Yang B. Review of optimization methodologies in global and China trajectory optimization competitions.
  Progress in Aerospace Sciences 102 (2018) 60–75. DOI 10.1016/j.paerosci.2018.07.004. https://www.sciencedirect.com/science/article/abs/pii/S0376042118300435
* [ZZ25] Zhang Z., Guo X., Wu D., Baoyin H., Li J., Topputo F. Global Optimality in Multi-Flyby Asteroid Trajectory Optimization: Theory and
  Application Techniques. arXiv 2508.02904 (2025); J. Guidance, Control, and Dynamics, DOI 10.2514/1.G009335. https://arxiv.org/abs/2508.02904
* [GT4] GTOC Portal, GTOC 4 – Asteroids billiard: https://sophia.estec.esa.int/gtoc_portal/?page_id=23 (problem description
  http://sophia.estec.esa.int/gtoc_portal/wp-content/uploads/2012/11/gtoc4_problem_description.pdf ; rankings
  http://sophia.estec.esa.int/gtoc_portal/wp-content/uploads/2012/11/ACT-RPT-MAD-GTOC4-ranks.pdf ; results summary
  http://sophia.estec.esa.int/gtoc_portal/wp-content/uploads/2014/11/gtoc4_summary_of_results.pdf)
* [GT4D] GTOC4 asteroid data: http://sophia.estec.esa.int/gtoc_portal/wp-content/uploads/2012/11/gtoc4_problem_data.txt (1435 bodies)
* [ALT] I. Althöfer, Combinatorial Space Trajectories to Comets and Asteroids (GTOC4 post-competition records): https://althofer.de/space-trajectories.html
* [ATK] ATK (Aerospace Tool Kit, NUDT): https://www.osredm.com/atknudt/atk/about
