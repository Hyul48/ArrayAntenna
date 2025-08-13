# ArrayAntenna

배열 안테나 실험 레포 (Julia). ULA 스티어링, Dolph–Chebyshev 가중치, 배열인자(Array Factor) 시뮬.

## Features
- Uniform/Binomial/Dolph–Chebyshev 가중치
- 조향(steering) & 빔폭/사이드로브 조절
- 패턴/AF 플로팅



## 📐 배열 가중치: Uniform / Binomial / Dolph–Chebyshev (ULA 기준)

### ULA 기본 정의
- 요소 수: $N$, 간격: $d$, 파수: $k=\frac{2\pi}{\lambda}$
- 조향각: $\theta_0$ (broadside면 $\theta_0=0$)
- 진행 위상(조향용): $\beta_n=-\,n\,k d \sin\theta_0$
- 편의 변수:
  
  $$
  \psi(\theta)=k d\big(\sin\theta-\sin\theta_0\big)
  $$
  
- 배열 인자(Array Factor):
  
  $$
  AF(\theta)=\sum_{n=0}^{N-1} w_n\,e^{\,j\,n\,\psi(\theta)}
  $$
  
  (중심대칭 배열이면 $w_n=w_{N-1-n}$.)


---

### 1) Uniform 가중치
- 정의: $w_n=1\ (n=0,\dots,N-1)$
- 닫힌꼴(Dirichlet kernel):
  
  $$
  AF_{\text{uni}}(\theta)
  = e^{\,j\frac{(N-1)\psi}{2}}\,
    \frac{\sin\!\big(\tfrac{N\psi}{2}\big)}{\sin\!\big(\tfrac{\psi}{2}\big)}
  $$
  
- 특징: 구현이 가장 단순, **사이드로브 $\approx -13.26$ dB**, 메인로브가 가장 좁음.  
  $d>\lambda/2$면 그레이팅 로브 주의.


---

### 2) Binomial(바이노미얼) 가중치
- 정의: $w_n=\binom{N-1}{n}$ (파스칼 삼각형 계수)
- AF 닫힌꼴:
  
  $$
  \begin{aligned}
  AF_{\text{bin}}(\theta)
  &= \sum_{n=0}^{N-1}\binom{N-1}{n}e^{\,j n\psi} \\
  &= \big(1+e^{\,j\psi}\big)^{N-1} \\
  &= e^{\,j\frac{(N-1)\psi}{2}}\;\Big[2\cos\!\big(\tfrac{\psi}{2}\big)\Big]^{N-1}
  \end{aligned}
  $$
  
- 특징: **이론상 사이드로브 0**(broadside). 대신 메인로브가 넓음(해상도 저하).


---

### 3) Dolph–Chebyshev(돌프–쳬비셰프) 가중치
- 목표: **사이드로브를 모두 동일 레벨**(equiripple)로 맞추면서 주엽 폭 최소화.
- Chebyshev 다항식(1종):
  
  $$
  T_m(x)=
  \begin{cases}
  \cos\!\big(m\arccos x\big), & |x|\le 1,\\[4pt]
  \cosh\!\big(m\,\operatorname{arcosh}x\big), & x\ge 1.
  \end{cases}
  $$
  
- 원하는 사이드로브 레벨(선형비) $1/R$을 정하면
  
  $$
  R=10^{\frac{\text{SLL(dB)}}{20}},\qquad
  x_0=\cosh\!\Big(\frac{1}{N-1}\operatorname{arcosh}R\Big)\ (\ge 1)
  $$
  
- 정규화된 AF 크기(피크=1, 사이드로브=$1/R$):
  
  $$
  \boxed{
  \big|AF_{\text{cheb}}(\theta)\big|
  = \frac{\big|T_{\,N-1}\!\big(x_0\,\cos(\psi(\theta)/2)\big)\big|}
         {T_{\,N-1}(x_0)}
  }
  $$

- 계수 비교를 통해서 가중치를 얻을 수 있음 $|A_1T_1(x) + A_2T_2(x)...| = |T_n(ax + b)|$

- 실제 $w_n$ 계산은 위 식의 $\cos(m\psi)$ 전개(또는 DFT/IDFT)로 계수 추출해서 대칭 $(w_n=w_{N-1-n})$을 맞추어 얻는다.  
  $R$을 키우면(SLL 더 낮춤) 사이드로브는 내려가고 메인로브는 넓어진다.


---

### 조향(steering)
진폭 가중 $w_n$은 그대로 두고 위상을

$$
\beta_n=-\,n\,k d\sin\theta_0
$$

로 주면 메인빔이 $\theta_0$로 향한다. 스캔 범위를 넓히려면 보통 $d\le \lambda/2$.


---

### 🔎 한눈 비교

| 가중치 | $w_n$ | 핵심 $AF(\theta)$ 식(요지) | 사이드로브 | 메인로브 폭 |
|---|---|---|---|---|
| Uniform | $1$ | $e^{j\frac{(N-1)\psi}{2}}\dfrac{\sin(N\psi/2)}{\sin(\psi/2)}$ | $\approx -13.26$ dB | 가장 좁음 |
| Binomial | $\binom{N-1}{n}$ | $e^{j\frac{(N-1)\psi}{2}}\,[2\cos(\psi/2)]^{N-1}$ | 0(이론상) | 가장 넓음 |
| Dolph–Chebyshev | (전개로 산출) | $\dfrac{\big\|T_{N-1}\!\big(x_0\cos(\psi/2)\big)\big\|}{T_{N-1}(x_0)}$ | 설계값 $1/R$ | SLL 고정 시 최소 |

> 기호: $\psi(\theta)=k d(\sin\theta-\sin\theta_0)$, $k=2\pi/\lambda$, $x_0=\cosh\!\big(\tfrac{1}{N-1}\operatorname{arcosh}R\big)$.

---

### 렌더링 체크리스트
- 이 섹션을 README.md에 붙일 때 **코드블록(\`\`\`)은 제외**해야 GitHub에서 수식이 렌더링됨.
- GitHub는 인라인 $…$, 블록 $$…$$ 수식을 지원. 붙여 넣고 **리포 페이지 새로고침**.
