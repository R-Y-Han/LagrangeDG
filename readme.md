# A Lagrangian DG hyperdynamic method

***Xiaodong Liu, Nathaniel R. Morgan, Donald E. Burton***

***

## ���Ʒ���

$$
\begin{equation}
\begin{split}
\frac{d\rho}{dt} & = - \rho \nabla \cdot \textbf{u}\\
\rho \frac{d\textbf{u}}{dt} & = \nabla \cdot \sigma\\
\rho \frac{d\tau}{dt} & = \nabla \cdot (\sigma \cdot \textbf{u})
\end{split}
\end{equation}
$$
����$\tau = e + k$Ϊ�������ܣ�$\sigma$ΪӦ������������������Ϊ$\sigma = -p \textbf{I}$��

***

## ��ɢ

$$
\textbf{x} = \sum_p b_p \textbf{x}_p
$$
����$\textbf{x}_p$Ϊ��Ԫ�Ľڵ����꣬$b_p$Ϊ�κ�����

$$
\mathbb{U}(\mathbf{\xi},t) = \mathbf{\psi}(\xi) \cdot \mathbb{U}^k(t)\\
\mathbb{U} = [\mathbb{U}_{cm}, \frac{\partial \mathbb{U}}{\partial \xi}_{cm}, \frac{\partial \mathbb{U}}{\partial \eta}_{cm}]^T\\
\psi = [1, \xi - \xi_{cm}, \eta - \eta_{cm}]^T
$$

### �����غ�

$$
\rho(\xi,t) |J| (\xi,t) = \rho^0(\xi) |J|(\xi,t^0)
$$

> *�����ʼ�ܶ��м����Ҫ���⴦�����ֵ���ܻ����ڼ�ϵĵط�����ǿ�����غ��������⡣�����ʼ�ܶ�û�м�ϲ���Ҫ�ر���*

### ʱ�����

$$
(\mathbb{U}^k)^{s1} = (\mathbb{U}^k)^n + \Delta t M^{-1} \cdot \mathbb{R}^n\\
(\mathbb{U}^k)^{n+1} = \frac{1}{2}(\mathbb{U}^k)^n + \frac{1}{2}(\mathbb{U}^k)^{s1} + \frac{1}{2}\Delta t M^{-1} \cdot \mathbb{R}^{s1}
$$

$$
(\textbf{x}_p)^{s1} = (\textbf{x}_p)^n + \Delta t \cdot (\textbf{u}_p^*)^n\\
(\textbf{x}_p)^{n+1} = \frac{1}{2}(\textbf{x}_p)^n + \frac{1}{2}(\textbf{x}_p)^{s1} + \frac{1}{2}\Delta t \cdot (\textbf{u}_p^*)^{s1}
$$

### �����غ���ɢ

$$
\begin{split}
& \left[
\begin{array}{ccc}
m & 0 & 0\\
0 &  \sum_{g\in \Omega} (\rho \psi_2 \psi_2)_g |J|_g \Omega_g &  \sum_{g\in \Omega} (\rho \psi_2 \psi_3)_g |J|_g \Omega_g\\
0 & \sum_{g\in \Omega} (\rho \psi_3 \psi_2)_g |J|_g \Omega_g &  \sum_{g\in \Omega} (\rho \psi_3 \psi_3)_g |J|_g \Omega_g
\end{array}
\right] \cdot \frac{\Delta}{\Delta t}
\left[
\begin{array}{c}
\textbf{u}_{cm}\\
\frac{\partial \textbf{u}}{\partial \xi}_{cm}\\
\frac{\partial \textbf{u}}{\partial \eta}_{cm}
\end{array}
\right]\\
& = \left[
\begin{array}{c}
\sum_{i\in \omega(t)} a_i \textbf{n}_i \cdot \sigma^*_c\\
\sum_{i\in \omega(t)} \psi_{2_p} a_i \textbf{n}_i \cdot \sigma^*_c\\
\sum_{i\in \omega(t)} \psi_{3_p} a_i \textbf{n}_i \cdot \sigma^*_c
\end{array}
\right] - \left[
\begin{array}{c}
0\\
\sum_{g \in \Omega} (\nabla_\xi \psi_2 \cdot |J| J^{-1} \cdot \sigma )_g \Omega_g\\
\sum_{g \in \Omega} (\nabla_\xi \psi_3 \cdot |J| J^{-1} \cdot \sigma )_g \Omega_g
\end{array}
\right]
\end{split}
$$

### �����غ���ɢ

$$
\begin{split}
& \left[
\begin{array}{ccc}
m & 0 & 0\\
0 &  \sum_{g\in \Omega} (\rho \psi_2 \psi_2)_g |J|_g \Omega_g &  \sum_{g\in \Omega} (\rho \psi_2 \psi_3)_g |J|_g \Omega_g\\
0 & \sum_{g\in \Omega} (\rho \psi_3 \psi_2)_g |J|_g \Omega_g &  \sum_{g\in \Omega} (\rho \psi_3 \psi_3)_g |J|_g \Omega_g
\end{array}
\right] \cdot \frac{\Delta}{\Delta t}
\left[
\begin{array}{c}
\tau_{cm}\\
\frac{\partial \tau}{\partial \xi}_{cm}\\
\frac{\partial \tau}{\partial \eta}_{cm}
\end{array}
\right]\\
& = \left[
\begin{array}{c}
\sum_{i\in \omega(t)} a_i \textbf{n}_i \cdot \sigma^*_c \cdot \textbf{u}_p^*\\
\sum_{i\in \omega(t)} \psi_{2_p} a_i \textbf{n}_i \cdot \sigma^*_c \cdot \textbf{u}_p^*\\
\sum_{i\in \omega(t)} \psi_{3_p} a_i \textbf{n}_i \cdot \sigma^*_c \cdot \textbf{u}_p^*
\end{array}
\right] - \left[
\begin{array}{c}
0\\
\sum_{g \in \Omega} (\nabla_\xi \psi_2 \cdot |J| J^{-1} \cdot \sigma \cdot \textbf{u})_g \Omega_g\\
\sum_{g \in \Omega} (\nabla_\xi \psi_3 \cdot |J| J^{-1} \cdot \sigma \cdot \textbf{u})_g \Omega_g
\end{array}
\right]
\end{split}
$$

***

## Riemann����

$$
\textbf{F}_i^* =  a_i \textbf{n}_i \cdot \sigma_c^* = a_i \textbf{n}_i \cdot \sigma_c + \mu_c |\textbf{n}_i \cdot \textbf{e}_c| a_i (\textbf{u}_p^* - \textbf{u}_c)
$$

����$\textbf{e}_c = \frac{\textbf{u}_p^*- \textbf{u}_c}{\Vert \textbf{u}_p^* - \textbf{u}_c \rVert}$����Ϊ$\textbf{u}_p^*$��ʱ��֪����������$\hat{\textbf{e}}_c = \frac{\textbf{u}_p^{avg}- \textbf{u}_c}{\lVert \textbf{u}_p^{avg} - \textbf{u}_c \rVert}$���棬$\textbf{u}_p^{avg}$Ϊ$\textbf{u}_c$��ƽ����

$$
\textbf{u}_p^* = \frac{\sum_{i \in p} \mu_c |\textbf{n}_i \cdot \hat{\textbf{e}}_c| a_i \textbf{u}_c - a_i \textbf{n}_i \sigma_c}{\sum_{i \in p}\mu_c |\textbf{n}_i \cdot \hat{\textbf{e}}_c|a_i}
$$

> *��Ϊ���㶯���غ�ĳ��������$\sum_i \textbf{F}^*_i = 0$������$\textbf{F}^*_i$�����е�$\textbf{e}_c$ҲҪ��$\hat{\textbf{e}}_c$����*
> $\sigma_c$�ڽڵ㴦ȡֵ��$\mu_c$�ڵ�Ԫ����ȡֵ

***

## ������

$\mathbb{U}^{lim}(\xi_p) = \overline{\mathbb{U}} + \phi (\mathbb{U}(\xi_p) - \overline{U})$

����
$\phi = 
\begin{cases}
\min(1,\frac{\alpha (\overline{\mathbb{U}}^{max}-\overline{\mathbb{U}})}{\mathbb{U}(\xi_p)-\overline{\mathbb{U}}}) \quad \mathbb{U}(\xi_p) > \overline{\mathbb{U}}\\
\min(1,\frac{\alpha (\overline{\mathbb{U}}^{min}-\overline{\mathbb{U}})}{\mathbb{U}(\xi_p)-\overline{\mathbb{U}}}) \quad \mathbb{U}(\xi_p) < \overline{\mathbb{U}}\\
1 \quad else
\end{cases}
$

***

## �㷨

**���㷨��û�����������*

> **t=0ʱ**
>
>> + ��������ͳ�ֵ
>> + ���㵥Ԫ���ĺ���������
>
> **$t_n$ʱ**
>
>> + ����ڵ��ٶ�$\textbf{u}_p^*$ ������$\mathbb{U}^n$��$\textbf{x}_p^n$)
>> + ������������$(\textbf{x}_p)^{s1} = (\textbf{x}_p)^n + \Delta t \cdot (\textbf{u}_p^*)^n$��ͬʱ��������������
>> + RK��һ����$(\mathbb{U}^k)^{s1} = (\mathbb{U}^k)^n + \Delta t M^{-1} \cdot \mathbb{R}^n$
>> + ����$(\textbf{u}_p^*)^{s1}$ ������$\mathbb{U}^{s1}$��$\textbf{x}_p^{s1}$��
>> + ��������ڵ�����$(\textbf{x}_p)^{n+1} = \frac{1}{2}(\textbf{x}_p)^n + \frac{1}{2}(\textbf{x}_p)^{s1} + \frac{1}{2}\Delta t \cdot (\textbf{u}_p^*)^{s1}$������������������
>> + RK�ڶ���$(\mathbb{U}^k)^{n+1} = \frac{1}{2}(\mathbb{U}^k)^n + \frac{1}{2}(\mathbb{U}^k)^{s1} + \frac{1}{2}\Delta t M^{-1} \cdot \mathbb{R}^{s1}$

***
***
***

## Test problems

***

### Shockless Noh

#### configuration

> $\gamma = 5/3$

Initial condition:

> $\rho^0 = 1$
>
> $u_x^0 = - x^0$
>
> $u_y^0 = - y^0$
>
> $e^0 = 1$

Boundary condition:
> &emsp; velocity of boundary nodes at any time are fixed:
>
> $u_x(x_b) = - x_b^0$, $u_y(y_b) = - y_b^0$

Analytic solution:

> $\rho = \rho^0(1-t)^q$
>
> $u = u^0$
>
> $v = v^0$
>
> $e = e^0(1-t)^{-q(\gamma - 1)}$
>
> $q = 2$

#### numerical result

��$t=0.6$ʱ�����ܵ����Ϊ��

|mesh          | $L^2$ error | order |
|:--:          |:-----------:|:-----:|
| 5$\times$5   |   4.545e-3  |   -   |
| 10$\times$10 |   1.133e-3  | 2.004 |
| 20$\times$20 |   2.830e-4  | 2.001 |
| 40$\times$40 |   7.072e-5  | 2.001 |
| 80$\times$80 |   1.768e-5  | 2.000 |

ѹ����ͼƬ��

![pressure](shockless_Noh/output/p.png)

***

### Taylor-Green vortex

#### configuration

> $\gamma = 7/5$

Initial condition:

> $\rho^0 = 1$
>
> $u_x^0 = \sin(\pi x) \cos(\pi y)$
>
> $u_y^0 = - \cos(\pi x) \sin(\pi y)$
>
> $p^0 = \frac{1}{4} [ \cos(2\pi x) + \cos(2\pi y)] + 1$
>
> with a source term in energy equation:
>
> $e_{s} = \frac{\pi}{4 ( \gamma - 1)}[\cos(3\pi x)\cos(\pi y) - \cos(3\pi y) \cos(\pi x)]$

Boundary condition:

> $u_n = 0$

Analytic solution:

> steady

#### numerical results

$t= 0.2$ʱѹ�����Ϊ��

|mesh          | $L^2$ error | order |
|:--:          |:-----------:|:-----:|
| 10$\times$10 |   4.286e-2  |   -   |
| 20$\times$20 |   1.265e-2  | 1.760 |
| 40$\times$40 |   3.589e-3  | 1.817 |
| 80$\times$80 |   1.146e-3  | 1.647 |

$t=0.2$ʱ��$40\times40$����ѹ�����£�

![t=0.2pressure](Taylor-Green_vortex/output/t=0.2_40.png)

$t=0.75$ʱ��$40\times40$�������£�

![t=0.75](Taylor-Green_vortex/output/t=0.75dirichlet.png)

֮ǰcorner force������Ϊ0ʱ�ᵼ�����´���

![t=0.75error](Taylor-Green_vortex/output/figure%201.%20Taylor-Green%20at%20t=0.75.png)

![t=0.75error2](Taylor-Green_vortex/output/figure%202.%20error%20near%20boundary%20.png)

***

### 1D Sod shock tube

#### configuration

> $\gamma = 5/3$

Initial conditions:

> $(\rho^0, u_x^0, u_y^0, p^0)_L = (1,0,0,1)$
>
> $(\rho^0, u_x^0, u_y^0, p^0)_R = (0.125,0,0,0.1)$

Boundary conditions:

> $u_n = 0$

#### numerical result

> ���������������Maire�Ľڵ�ⷨ����
> $\textbf{F}^*_{i,maire} = a_i \textbf{n}_i \sigma^*_c = a_i \textbf{n}_i \sigma_c +a_i \mu_c (\textbf{n}_i \otimes \textbf{n}_i) (\textbf{u}_p^* - \textbf{u}_c)$
> ��Burton���������ڣ�
> $$\textbf{F}^*_{i,burton} = a_i \textbf{n}_i \sigma^*_c = a_i \textbf{n}_i \sigma_c +a_i \mu_c (\textbf{e}_c \otimes \textbf{e}_c) (\textbf{u}_p^* - \textbf{u}_c)\\
\quad \quad \quad = a_i \textbf{n}_i \sigma_c +a_i \mu_c  \textbf{e}_c |\textbf{u}_p^* - \textbf{u}_c|$$
> Burton�ⷨ����$\textbf{e}_c$�а�����$\textbf{u}_p^*$�������ٶȳ�ֵ������0ʱΪ�õ��ٶ���Ҫ��$\textbf{e}_c$���޸ģ���ͬ�޸ĵķ��������һ��

��$t=0.2$ʱ������£�

![pressure](sod_shock_tube/output/p,maire,cellcentered_rho.png)

![internal_energy](sod_shock_tube/output/e,maire,cellcentered_rho.png)

![density](sod_shock_tube/output/rho,maire,cellcentered_rho.png)

![velocity](sod_shock_tube/output/velocity,maire,cellcentered_rho.png)

***

### Sedov blast

#### configuration

>$\gamma = 7/5$

Initial conditions:

> $(\rho^0, u_x^0, u_y^0, p^0) = (1,0,0,10^{-6})$
>
> energy source at the origin:
>
> $p_O = (\gamma - 1)\rho_O \frac{E_O}{V_O}$, $E_O$Ϊ�ͷŵ������ܣ�$V_O$Ϊ���ļ�������������������Ϊ0.5

#### numerical result

��$t=1.5$ʱѹ����ͼ��Ϊ��

![pressuret=1.5](sedov/output/t=1.5.png)

�����ܶ�ͼ��Ϊ��
![density]

### Noh

#### configuration

> $\gamma = 5/3$

Initial condition��

> ��ʼ�ٶ�Ϊָ�����ĵĵ�λ��������ʼ�ܶȺ�Ϊ1����ʼѹ��Ϊ1e-6

boundary condition:

> $p_b \equiv 10^{-6}$

#### numerical results

$t=0.6$ʱ����ѹ�����£�

![pressure](Noh/output/t=0.6.png)