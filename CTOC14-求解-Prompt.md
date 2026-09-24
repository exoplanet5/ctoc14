你是一名航天器轨道设计、最优控制、组合优化和科学计算专家。请帮助我系统求解第十四届全国空间轨道设计竞赛 CTOC14 题目甲：

**“近地小行星防御多航天器遍历探测轨道设计与优化”**

我会提供两个附件：

1. `CTOC14_problem.pdf`
2. `MEA.txt`

请始终以这两个文件中的定义、参数、动力学、约束、评分函数和提交格式为最终权威，不要自行修改题意。

---

# 1. 最终目标

建立一个完整、可迭代、可验证的 CTOC14 求解框架。

核心任务是在 2030-01-01 起的 15 年绝对任务窗口内，用若干艘连续低推力航天器飞越探测 300 颗给定近地小行星，并最小化：

\[
J=\sum_{i=1}^N J_i+N_{\rm miss}
\]

其中

\[
J_i=1+x_i+x_i^2,
\qquad
x_i=\frac{m_{0,i}-600}{1400}.
\]

最终还需要考虑题目规定的提交时间系数，但轨迹优化阶段首先最小化原始任务代价 \(J\)。

我们的目标不是只做一个能够工作的 baseline，而是建立一个有潜力用于竞赛排名的完整 solver。

---

# 2. 必须严格遵守的题目模型

## 2.1 坐标系

全部状态均使用：

Heliocentric Ecliptic Inertial, HEI

- 原点：太阳质心
- xy 平面：J2000.0 黄道面
- x 轴：J2000.0 春分点
- z 轴：黄道北极方向

所有位置单位统一使用 km。

所有速度单位统一使用 km/s。

动力学积分内部时间统一使用 s。

---

# 2.2 航天器动力学

航天器满足日心二体 + 连续推力：

\[
\ddot{\mathbf r}
=
-\frac{\mu_\odot}{r^3}\mathbf r
+
\frac{\mathbf T}{m}.
\]

其中

\[
0\leq \|\mathbf T\|\leq0.5\ {\rm N}.
\]

质量变化：

\[
\dot m
=
-\frac{\|\mathbf T\|}{I_{sp}g_0}.
\]

参数：

\[
\mu_\odot
=
1.32712440018\times10^{11}\ {\rm km^3/s^2}
\]

\[
I_{sp}=4000\ {\rm s}
\]

\[
g_0=9.80665\ {\rm m/s^2}
\]

\[
T_{\max}=0.5\ {\rm N}.
\]

注意动力学中推力产生的加速度必须正确处理 N、kg、km/s² 的单位换算。

---

# 2.3 航天器质量

干质量：

\[
m_{\rm dry}=600\ {\rm kg}.
\]

初始质量：

\[
600\leq m_0\leq2000\ {\rm kg}.
\]

燃料质量：

\[
m_{\rm fuel}=m_0-600.
\]

任意时刻均必须满足：

\[
m(t)\geq600\ {\rm kg}.
\]

---

# 2.4 地球出发条件

任务起点：

2030-01-01 00:00:00 UTC

即：

\[
{\rm MJD}=62502.0.
\]

每艘航天器发射时间可独立选择：

\[
t_{\rm launch}\geq0.
\]

且必须在整个 15 年窗口以内。

发射时：

\[
\mathbf r_{\rm SC}
=
\mathbf r_\oplus.
\]

速度满足：

\[
\left\|
\mathbf v_{\rm SC}
-
\mathbf v_\oplus
\right\|
\leq4\ {\rm km/s}.
\]

这个 \(v_\infty\) 不进入燃料代价，因此应视为一个极其重要的“免费初始速度控制量”。

方向可任意选择。

需要在算法中充分利用这一点。

---

# 2.5 地球和小行星星历

所有天体只按太阳中心二体问题传播。

地球轨道根数：

历元：

MJD 60676.0

轨道根数直接从题目 PDF 读取。

300 颗小行星的轨道根数由：

`MEA.txt`

提供。

其统一历元：

\[
{\rm MJD}=61200.0.
\]

必须自己实现经典轨道根数：

\[
(a,e,i,\Omega,\omega,M_0)
\]

到任意时刻 Cartesian state：

\[
(\mathbf r,\mathbf v)
\]

的二体传播。

需要高精度解决 Kepler 方程。

不要使用外部在线星历替代题目星历。

---

# 2.6 小行星探测条件

只要求飞越。

当：

\[
\left\|
\mathbf r_{\rm SC}(t)
-
\mathbf r_j(t)
\right\|
\leq1000\ {\rm km}
\]

即认为完成目标 \(j\) 的探测。

非常重要：

不要求：

\[
\mathbf v_{\rm SC}
=
\mathbf v_j.
\]

因此这里是 flyby，不是 rendezvous。

小行星不提供任何引力辅助。

飞越前后航天器状态连续。

重复访问同一个目标不重复计分。

---

# 2.7 任务时间

所有事件都必须满足：

\[
0\leq t\leq5478.75\ {\rm days}.
\]

即：

15 年。

超出窗口的飞越不计覆盖。

---

# 3. 评分函数的算法含义

每艘航天器：

\[
J_i
=
1+x_i+x_i^2.
\]

因此：

- 发射一艘干质量航天器：成本 1
- 满燃料 2000 kg：成本 3
- 每漏掉一个 asteroid：成本同样为 1

所以 solver 的本质不是单纯最小燃料，而是：

**最大化每艘航天器带来的有效目标覆盖价值。**

对于一条路线 \(R\)，如果覆盖 \(K_R\) 个新目标，估计燃料归一化量为 \(x_R\)，可以定义 route utility：

\[
U(R)
=
K_R
-
(1+x_R+x_R^2).
\]

这是路线生成和筛选的重要评价指标。

不要人为规定航天器数量 \(N\)。

让组合优化最终自动决定：

- 发多少艘航天器
- 每艘访问哪些 asteroid
- 哪些极难目标可以战略性放弃

---

# 4. 总体算法架构

不要直接把：

- 300 个 asteroid assignment
- 所有发射时间
- 所有 encounter time
- 所有低推力控制

一次性放进一个巨大 NLP。

这样规模太大且高度非凸。

采用分层求解：

\[
\boxed{
\text{Ephemeris}
\rightarrow
\text{Ballistic Atlas}
\rightarrow
\text{State-expanded Graph}
\rightarrow
\text{Route Search}
\rightarrow
\text{Set Cover}
\rightarrow
\text{Low-thrust Refinement}
\rightarrow
\text{Submission Validator}
}
\]

---

# 5. Phase 0：基础动力学与验证

首先实现：

```text
ctoc14/
    constants.py
    ephemeris.py
    kepler.py
    lambert.py
    dynamics.py
    propulsion.py
    data.py
    tests/
```

实现：

```python
state_from_elements(elements, mjd)
earth_state(mjd)
asteroid_state(id, mjd)
```

以及：

```python
propagate_two_body(...)
```

要求：

- 使用 double precision
- 单位体系清晰统一
- 对 300 个目标批量传播足够快
- 尽可能 vectorize

---

# 6. 第一项严格验证

利用题目 PDF 的 SC1 示例：

2030-01-01 从地球发射，之后弹道飞向 asteroid 174。

反解 Earth → 174 Lambert。

验证应能够重现题目示例的大致结果：

\[
v_\infty\simeq4\ {\rm km/s}.
\]

同时验证：

- Earth state
- asteroid 174 state
- Lambert solver
- 时间转换
- HEI 坐标方向
- 单位体系

这必须成为 automated regression test。

如果无法复现示例，不要继续后续优化。

---

# 7. Phase 1：Earth → asteroid ballistic atlas

这一阶段完全不使用电推进。

目的：

扫描所有可能的：

\[
Earth(t_L)\rightarrow Asteroid_j(t_A)
\]

Lambert 轨迹。

变量：

\[
t_L,\quad t_A.
\]

约束：

\[
t_A>t_L
\]

以及：

\[
v_\infty
=
\left\|
\mathbf v_{\rm Lambert}(t_L)
-
\mathbf v_\oplus(t_L)
\right\|
\leq4\ {\rm km/s}.
\]

需要同时考虑 Lambert 的合理分支：

- short-way
- long-way
- 必要时 multi-revolution

但第一版可以从 zero-revolution 开始。

---

# 8. Ballistic atlas 的搜索方式

第一轮使用粗网格：

发射时间间隔例如：

30–90 days。

飞行时间：

30 days 到若干年。

第一版可限制在：

50–1500 days。

对每个 asteroid 搜索。

之后对 promising windows 用局部优化精化：

\[
\min_{t_L,t_A} v_\infty.
\]

每颗 asteroid 不只保存一个窗口。

保存前：

50–200 个

distinct arrival states。

定义：

```python
class ArrivalState:
    asteroid_id
    launch_time
    arrival_time
    r_arr
    v_arr_sc
    v_arr_ast
    vinf
```

注意最关键的状态变量是：

\[
\mathbf v_{\rm SC}(t_A)
\]

而不是 asteroid 本身的速度。

---

# 9. 为什么必须保存 state-expanded nodes

对于同一目标 \(A\)，可能存在很多：

\[
(A,t_A,\mathbf v_{\rm SC})
\]

状态。

虽然都飞越 A，但之后去 B 的难度可能完全不同。

所以后续图不是：

\[
A\rightarrow B
\]

而是：

\[
(A,t,\mathbf v)
\rightarrow
(B,t',\mathbf v').
\]

这叫做：

state-expanded temporal graph。

---

# 10. Phase 2：asteroid-to-asteroid impulsive surrogate

对已有状态：

\[
S_A=
(A,t_A,\mathbf r_A,\mathbf v_{\rm SC,A})
\]

枚举未来 target B 和：

\[
t_B>t_A.
\]

求 Lambert：

\[
\mathbf r_A(t_A)
\rightarrow
\mathbf r_B(t_B).
\]

Lambert 所需 departure velocity：

\[
\mathbf v_L.
\]

定义 surrogate：

\[
\Delta v_{\rm proxy}
=
\|
\mathbf v_L-\mathbf v_{\rm SC,A}
\|.
\]

Lambert 到达速度：

\[
\mathbf v_{\rm arr,L}
\]

作为下一状态的：

\[
\mathbf v_{\rm SC,B}.
\]

这样形成：

\[
S_A\rightarrow S_B.
\]

注意：

这个 \(\Delta v_{\rm proxy}\) 不是真正连续低推力最优值。

它只是一个 cheap heuristic，用来评价边是否值得深入优化。

---

# 11. Edge pruning

绝不能保存完整 dense graph。

对于每个 state，只保存最有希望的若干 successor。

可以结合：

\[
\Delta v_{\rm proxy}
\]

flight time

inclination / geometry

remaining mission time

future reachable target density

建立排序。

例如定义：

\[
C_{AB}
=
w_1\Delta v_{\rm proxy}
+
w_2\Delta t
-
w_3N_{\rm future}.
\]

每个节点只保留：

20–100

条最佳 edge。

---

# 12. 连续推力粗略可实现性模型

利用：

\[
a_T=\frac{T_{\max}}{m}
\]

估计实现 impulsive \(\Delta v\) 所需推力时间。

粗略：

\[
t_{\rm thrust}
\approx
\frac{\Delta v_{\rm proxy}}{T/m}.
\]

进一步可用 rocket equation 估算 fuel：

\[
m_f
=
m_i
\exp
\left(
-\frac{\Delta v}{I_{sp}g_0}
\right).
\]

这只是 surrogate。

用来：

- prune
- route ranking
- fuel estimate

最终必须由真正低推力 solver 重新优化。

---

# 13. Phase 3：route generation

现在需要搜索：

\[
Earth
\rightarrow
A_1
\rightarrow
A_2
\rightarrow
...
\rightarrow
A_K.
\]

建议实现至少以下一种：

1. Beam Search
2. GRASP
3. ALNS

第一版优先 Beam Search。

每个 beam state 保存：

```python
RouteState:
    time
    position
    velocity
    estimated_mass
    visited_mask
    asteroid_sequence
    encounter_times
    launch_time
    vinf_vector
    estimated_cost
```

---

# 14. Beam search 评分函数

不要只最小化累计 \(\Delta v\)。

建议类似：

\[
F
=
-N_{\rm visited}
+
\lambda_1 J_{\rm est}
+
\lambda_2\Delta v_{\rm remaining}
+
\lambda_3\frac{t}{T_{\max}}
-
\lambda_4N_{\rm future}.
\]

或者最大化：

\[
F
=
N_{\rm visited}
-
\lambda J_{\rm spacecraft}
+
\eta N_{\rm future}.
\]

需要实验不同权重。

核心目标：

找到高 target-density route。

---

# 15. 特别搜索 ballistic chains

建立专门算法搜索：

\[
Earth\rightarrow A\rightarrow B\rightarrow C
\]

几乎不使用低推力的路线。

因为：

\[
v_\infty\leq4\ {\rm km/s}
\]

是免费控制量。

如果一艘 600 kg 航天器不使用任何推进，却能覆盖：

2 个及以上目标，

这通常是极高价值路线。

因此需要实现：

```python
search_ballistic_chains()
```

允许对 Earth departure：

\[
t_L,\quad
v_\infty,\quad
RA,\quad Dec
\]

或者直接 Lambert-derived departure states

进行大规模搜索。

---

# 16. 可选的 Monte Carlo ballistic screening

另外建立一个随机搜索模式：

随机生成：

\[
t_L
\]

和

\[
\mathbf v_\infty,
\qquad
\|\mathbf v_\infty\|\leq4.
\]

把航天器作为纯日心二体轨道传播 15 年。

与此同时预计算所有 asteroid：

\[
\mathbf r_j(t_k).
\]

使用：

- KD-tree
- Ball tree
- spatial hashing

查询轨迹附近 asteroid。

筛出 near miss，例如：

\[
d<10^5\ {\rm km}
\]

甚至：

\[
10^6\ {\rm km}.
\]

再对这些 candidate 做局部优化，把距离压到：

\[
1000\ {\rm km}.
\]

---

# 17. Phase 4：route pool

大量搜索之后建立 route database：

```python
Route:
    id
    launch_time
    vinf
    asteroid_sequence
    encounter_times
    estimated_fuel
    estimated_initial_mass
    estimated_cost
    asteroid_mask
```

目标生成：

至少数千条，

理想情况下：

\[
10^4-10^6
\]

条 candidate route。

但需要 dominance pruning。

---

# 18. Route dominance

如果两个 route：

\[
R_1,\ R_2
\]

满足：

\[
S_{R_1}
\supseteq
S_{R_2}
\]

且：

\[
c_{R_1}\leq c_{R_2},
\]

则 \(R_2\) 被支配，可以删除。

类似地，对相同 asteroid set：

只保留若干：

- 最低 fuel
- 最早 finish
- 最有后续扩展潜力

Pareto optimal routes。

---

# 19. Phase 5：global route selection

建立 weighted set-cover / prize-collecting set-cover。

变量：

\[
y_r\in\{0,1\}
\]

表示是否选路线 r。

\[
z_j\in\{0,1\}
\]

表示 asteroid j 是否遗漏。

优化：

\[
\min
\sum_r
c_r y_r
+
\sum_j z_j.
\]

约束：

\[
\sum_{r:j\in S_r} y_r
+
z_j
\geq1.
\]

可使用：

- scipy milp
- OR-Tools
- HiGHS
- Gurobi（如环境存在）

优先让代码不依赖商业 solver。

---

# 20. 不要求避免 route overlap

多个航天器可以重复飞越相同 asteroid。

只是重复访问不额外计分。

因此 set cover 不需要 exact cover。

允许：

\[
\sum_r y_r > 1
\]

覆盖同一 target。

---

# 21. Hard targets 单独分析

在 route generation 前先分析 300 个目标的轨道统计。

重点列出：

- high inclination
- retrograde
- large semimajor axis
- 2030–2045 中始终远离内太阳系
- 很难通过 Earth ballistic injection 到达的对象

尤其单独分析：

- asteroid 131
- asteroid 144

不要让极少数异常目标严重干扰主体 298 颗目标的 graph search。

允许建立：

```text
main_target_pool
hard_target_pool
```

对 hard targets 设计专用任务。

---

# 22. asteroid 131 的策略

不要因为它是 retrograde 就错误地要求速度匹配。

仍然只是 position intercept。

因此研究：

- 2030 年尽早发射
- 大出平面轨迹
- 长时间连续推力
- 是否可在尚未远离太阳前拦截

求解变量：

\[
t_L
\]

\[
\mathbf v_\infty
\]

连续推力。

可以单独做一个：

Earth → 131

最小初始质量问题。

---

# 23. asteroid 144 的策略

这是超大半长轴、高偏心率目标。

先传播整个 2030–2045 窗口的：

\[
r(t).
\]

然后设计：

fast heliocentric escape trajectory。

同样不需要匹配其速度。

研究：

\[
Earth\rightarrow144
\]

是否能在 15 年内实现 position intercept。

如果可实现，求：

\[
m_{0,\min}.
\]

再根据：

\[
1+x+x^2
\]

判断覆盖它是否比：

\[
N_{\rm miss}=1
\]

更划算。

注意：

如果一艘专船成本 >1，却只覆盖一个 asteroid，从评分上通常不划算。

除非：

- 它还能串其它目标
- 或总成本增量因为某种设计小于漏掉目标的 penalty

所以 hard-target 任务必须做 marginal score 分析。

---

# 24. Phase 6：真正的 low-thrust refinement

当 asteroid sequence 已经固定：

\[
Earth
\rightarrow
A_1
\rightarrow
...
\rightarrow
A_K
\]

再做连续低推力 NLP。

优先实现：

direct multiple shooting

或者：

direct collocation。

第一版推荐：

direct multiple shooting。

---

# 25. Low-thrust 优化变量

至少包括：

\[
t_{\rm launch}
\]

\[
\mathbf v_\infty
\]

各 encounter time：

\[
t_1,\ldots,t_K
\]

以及 piecewise thrust control。

例如每段用：

\[
N_c
\]

个 control intervals。

每段控制：

\[
\mathbf T_k.
\]

约束：

\[
\|\mathbf T_k\|\leq0.5N.
\]

---

# 26. Encounter constraint

不要写 rendezvous：

错误：

\[
\mathbf r_{\rm SC}
=
\mathbf r_{\rm ast}
\]

同时

\[
\mathbf v_{\rm SC}
=
\mathbf v_{\rm ast}.
\]

正确：

\[
\|
\mathbf r_{\rm SC}(t_k)
-
\mathbf r_{A_k}(t_k)
\|
\leq1000\ {\rm km}.
\]

初期优化为了数值稳定可先设：

\[
\mathbf r_{\rm SC}
=
\mathbf r_{A_k}
\]

作为 equality。

得到可行解以后再允许：

1000 km

自由度。

---

# 27. Low-thrust objective

对固定路线，主要最小化：

\[
m_0
\]

或等价最小化 fuel。

因为：

\[
J_i
=
1+x+x^2
\]

对 \(m_0\) 单调增加。

最终使用实际：

\[
J_i
\]

作为评价。

---

# 28. 初始猜测

低推力 NLP 很依赖 initial guess。

用 Phase 2 的 Lambert impulsive route 生成 guess。

将 impulsive burn：

\[
\Delta \mathbf v
\]

替换成附近一段 finite burn。

例如按最大推力估计 burn duration：

\[
\tau
\sim
\frac{m\Delta v}{T}.
\]

方向初始设为：

\[
\hat{\mathbf T}
=
\frac{\Delta\mathbf v}{\|\Delta\mathbf v\|}.
\]

然后交给 nonlinear solver 自行优化。

---

# 29. 可使用的数值工具

Python 优先。

可以使用：

- numpy
- scipy
- numba
- joblib
- multiprocessing
- pandas
- matplotlib

优化器可选择：

- scipy.optimize
- scipy.optimize.least_squares
- scipy.optimize.minimize
- scipy.integrate.solve_ivp

如果安装了：

- casadi
- cyipopt

则优先考虑 CasADi + IPOPT 做连续推力 NLP。

但基础代码不要强依赖稀有商业软件。

---

# 30. 性能要求

大量 Lambert 搜索将成为瓶颈。

需要：

- vectorization
- parallel processing
- caching ephemerides
- precomputed asteroid state grids
- numba where useful

建立：

```python
EphemerisCache
LambertCache
```

避免重复传播。

---

# 31. 时间表示

代码内部建议统一：

```python
t_sec
```

表示：

从 MJD 62502.0 起的秒数。

需要函数：

```python
t_to_mjd(t_sec)
mjd_to_t(mjd)
```

以及：

```python
t_to_datetime(...)
```

所有 optimizer 内部尽量只使用秒或天的相对时间。

不要在核心动力学中混用：

MJD、day、second。

---

# 32. 最终 submission thrust representation

题目最后不是直接读取 optimizer 的 piecewise-constant control。

提交 Event=1 的推力样点后，组委会会使用：

3 阶滑动窗口 Lagrange interpolation

重新构造：

\[
\mathbf T(t).
\]

因此最终控制离散化必须显式适配题目这一规则。

需要单独实现：

```python
submission_thrust_interpolator()
```

完全复现官方规则。

---

# 33. 推力样点约束

同一 thrust arc 内：

Event=1 相邻采样时间必须满足：

\[
\Delta t\geq8640\ {\rm s}.
\]

并且不仅采样点：

\[
\|\mathbf T_k\|\leq0.5N
\]

还要保证 interpolation 后任意时刻：

\[
\|\mathbf T(t)\|\leq0.5N.
\]

三阶 Lagrange interpolation 可能 overshoot。

因此最终不要让节点直接打满：

\[
0.5N.
\]

需要保留 margin，例如先限制：

\[
0.47-0.49N
\]

或者对每个 interpolation interval 显式搜索最大：

\[
\|\mathbf T(t)\|.
\]

---

# 34. Submission validator

必须独立实现一个尽可能复刻官方的：

```python
validate_submission(path)
```

逐行检查：

- SC_ID
- Event sequence
- strictly increasing time
- launch condition
- \(v_\infty\)
- mass
- thrust
- thrust interpolation
- asteroid flyby
- mission window
- Event fields

并逐段从前一个 state 数值积分至下一个 row。

误差目标比官方更严格，例如：

position：

\[
<0.1\ {\rm km}
\]

velocity：

\[
<0.1\ {\rm m/s}
\]

mass：

\[
<0.001\ {\rm kg}.
\]

官方允许值再从 PDF 读取并作为 hard validation threshold。

---

# 35. 输出精度

最终文件至少：

```python
"%.17g"
```

或：

```python
"%.12e"
```

不要为了减小文件体积损失状态精度。

---

# 36. 软件目录结构

建议最终形成：

```text
ctoc14/
│
├── data/
│   └── MEA.txt
│
├── constants.py
├── time_utils.py
├── kepler.py
├── ephemeris.py
├── lambert.py
├── dynamics.py
├── propulsion.py
│
├── ballistic_atlas.py
├── asteroid_graph.py
├── beam_search.py
├── monte_carlo.py
├── route_pool.py
├── set_cover.py
│
├── low_thrust/
│   ├── transcription.py
│   ├── shooting.py
│   ├── controls.py
│   └── refine.py
│
├── submission/
│   ├── interpolate.py
│   ├── export.py
│   └── validator.py
│
├── analysis/
│   ├── target_statistics.py
│   ├── hard_targets.py
│   ├── route_statistics.py
│   └── score_analysis.py
│
├── cache/
├── results/
├── tests/
│
└── main.py
```

---

# 37. 第一阶段实际交付任务

不要一开始就写完所有模块。

首先完成以下 milestone。

## Milestone A

读取：

`MEA.txt`

输出 300 个 asteroid 的：

- a
- e
- i
- q
- Q
- orbital period

并绘制：

- a-e
- a-i
- q-i
- 2030–2045 heliocentric distance envelope

自动识别离群目标。

---

## Milestone B

实现题目专用 ephemeris。

验证 Earth 和 asteroid states。

特别利用 PDF sample：

Earth → asteroid 174

做 regression test。

---

## Milestone C

实现 Lambert solver 或采用经过验证的本地库。

搜索：

Earth → all 300 asteroids

ballistic flyby windows。

输出：

```text
asteroid_id
launch_time
arrival_time
TOF
vinf
arrival spacecraft velocity
```

每个 asteroid 保存多个不同窗口。

---

## Milestone D

统计：

有多少 asteroid 可以用：

\[
v_\infty\leq4{\rm km/s}
\]

纯弹道直接到达。

按：

- minimum vinf
- earliest reachable
- number of launch windows
- inclination
- semimajor axis

分析。

---

## Milestone E

从 ballistic arrival states 构造：

asteroid → asteroid Lambert transitions。

先只对：

10–50

个最容易目标实验。

实现 beam search，找：

2-target

3-target

5-target

甚至更长路线。

---

# 38. 每一步都必须给出可验证数字

不要只写框架。

每完成一个模块，都需要实际运行并输出：

- 测试结果
- numerical errors
- sample trajectories
- timing
- search statistics

例如：

```text
Earth->174 regression:
position error = ...
velocity error = ...
vinf = ...
```

以及：

```text
Ballistic atlas:
targets reachable = ...
Lambert solves = ...
CPU time = ...
```

---

# 39. 不要隐藏失败案例

任何：

- Lambert no solution
- optimizer failure
- collision with time boundary
- negative mass
- interpolation overshoot
- shooting divergence

都要记录。

使用 structured logging。

---

# 40. 优先研究的问题

当 baseline 建立后，依次研究：

### A

纯 ballistic route 能否自然覆盖多个目标？

### B

一次小幅速度修正是否能从某个 ballistic chain 切到另一个 chain？

### C

哪些 asteroid 是网络中的 high-degree hub？

### D

哪些 asteroid 是必须单独解决的 outlier？

### E

增加 100 kg 燃料平均能多覆盖多少目标？

### F

最优任务是否倾向：

- 很多廉价 spacecraft
- 中等数量、中等燃料
- 少量满燃料长链 spacecraft

### G

全覆盖是否真的值得？

因为漏掉：

\[
1\text{ target}
\]

只罚：

\[
1.
\]

对于昂贵的孤立目标，可能主动放弃才是最优策略。

---

# 41. 一个重要的 marginal-cost criterion

对于一艘新航天器：

\[
c=1+x+x^2.
\]

如果它只能新覆盖：

\[
K
\]

个 asteroid，

只有当：

\[
K>c
\]

时，相对于全部漏掉这些目标才严格降低 score。

在 route generation、hard target analysis 和 local improvement 中显式使用这一判断。

---

# 42. Local improvement

一旦有完整任务方案，再实现：

- route merge
- route split
- asteroid insertion
- asteroid removal
- route swap
- subsequence replacement
- launch-time shift
- encounter-time shift

例如：

如果：

SC1：

\[
A\rightarrow B\rightarrow C
\]

SC2：

\[
D\rightarrow E
\]

尝试：

\[
A\rightarrow B\rightarrow D\rightarrow E
\]

以及：

\[
C
\]

是否能插入其它 route。

这类似 vehicle-routing local search。

---

# 43. ALNS operators

后期可实现：

Destroy：

- random asteroid removal
- expensive asteroid removal
- route removal
- hard cluster removal

Repair：

- cheapest insertion
- regret-2 insertion
- regret-k insertion
- new spacecraft creation

评价始终使用真实：

\[
J
\]

或尽量可靠的 surrogate。

---

# 44. 允许不全覆盖

必须始终保留：

\[
N_{\rm miss}
\]

作为显式优化变量。

不要强制 300 全覆盖。

最终至少给出：

- best full-coverage solution
- best unconstrained-score solution

如果 full coverage 可行。

---

# 45. 最终需要维护 leaderboard

每次得到更好任务方案，记录：

```text
run_id
timestamp
J
N_spacecraft
N_covered
N_missed
total_fuel
max_fuel
mean_targets_per_SC
max_targets_per_SC
```

同时保存：

```text
solution.json
routes.csv
submission.txt
```

不要覆盖旧结果。

---

# 46. Reproducibility

所有 stochastic optimizer：

- Monte Carlo
- random restart
- GRASP
- ALNS

都显式记录 random seed。

配置全部存入：

```yaml
config.yaml
```

---

# 47. 可视化

至少输出：

1. 300 asteroid orbital elements 分布
2. target accessibility histogram
3. asteroid temporal graph
4. route length distribution
5. fuel vs asteroid count
6. spacecraft timeline
7. 单条 spacecraft 3D heliocentric trajectory
8. encounter chronology
9. overall coverage map
10. current best score evolution

---

# 48. 数值原则

任何时候都区分：

### Exact problem quantity

题目规定的真实动力学和 score。

### Surrogate

例如：

- Lambert \(\Delta v\)
- impulsive approximation
- simple thrust-time estimate

所有输出中清楚注明。

不要把 surrogate 误认为真实任务代价。

---

# 49. 最重要的建模认识

始终牢记：

这是：

**moving-target flyby routing**

而不是：

**multiple rendezvous mission**。

在 asteroid flyby：

\[
\mathbf v_{\rm SC}
\]

可以与 asteroid velocity 完全不同。

这一点是整个题目高效求解的核心。

---

# 50. 工作方式

请按“实现 → 运行 → 检验 → 分析 → 下一步”的方式工作。

不要一次输出几千行未经测试的代码。

每完成一个阶段：

1. 说明实现内容
2. 实际运行
3. 报告数字
4. 验证物理正确性
5. 保存结果
6. 再进入下一阶段

如果发现前面的算法假设不成立，应修改求解路线，而不是机械坚持初始方案。

---

# 51. 当前首先执行的任务

现在从 Phase 0 和 Phase 1 开始。

具体执行：

1. 读取 `CTOC14_problem.pdf`
2. 读取 `MEA.txt`
3. 自动提取并核对全部题目常数
4. 分析 300 asteroid 的轨道统计
5. 实现题目规定的二体星历
6. 用 asteroid 174 的 PDF 示例验证星历和 Lambert
7. 建立第一版 Earth → asteroid ballistic atlas
8. 输出所有 300 个目标的 minimum ballistic \(v_\infty\)
9. 统计 \(v_\infty\le4\) km/s 的目标
10. 识别 hard targets
11. 基于结果提出 Phase 2 graph-search 的具体参数

完成这些以后，再继续 asteroid-to-asteroid route search。

不要跳过验证直接进入低推力最优控制。