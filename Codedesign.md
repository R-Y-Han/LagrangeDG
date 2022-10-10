# Cell-centered Lagrange DG Code

***

## &nbsp; config.h

**����һЩ��������Ҫʹ�õĲ����ͺ�����**

+ ��������Ĳ���
  + ���񻮷ֵĴ�Сn, m
  + ��������x, y
  + �ռ䲽��$h_x, h_y$
  + ��ֹʱ��T
  + ����ʱ�䲽��dt
  + ����ı�����$\gamma$
  + EOS
+ �����趨
  + ���񡢽ڵ�ṹ�嶨��
    + ����������Ϣ�Ͳ���������
    + ��������Ͳο�������ͬһ������ṹ�壬���ο���Ԫ�ļ�����ϢҲ����������
+ ��ѧ��⹤��
  + [Gauss���](https://www.cnblogs.com/sunzhenwei/p/10848686.html)+Gauss�ڵ�
  + �κ���������ȷ���Ӳο���Ԫ�����������ӳ��
  + ����䳤�ȼ��㡢�ⷨ�������㣨ע�ⶨ��

**������Ϊ�ṹ��Ӧ������**

+ ��ţ�һά����ά����ע��ת����ʽ��
+ ���㣨��һ��˳��
+ ���ĵĲο�����
+ ������[���ṹ���ж��庯����](https://blog.csdn.net/weixin_45525272/article/details/106446669)
+ ���ڵ�Ԫ
+ Jacobi����
+ ��������
+ �����������꣨���ڼ���Jacobi����ΪJacobi����װ�ˣ�

**�ڵ���Ϊ�ṹ��Ӧ������**

+ ��ţ�һά����ά��ע��ת����ʽ��
+ ��������
+ �ڵ��ٶ�
+ ���ڵ�Ԫ
+ ���ڽڵ㣨����ص����ڼ���ڵ��ٶȵı߳����ⷨ��������
+ �Ƿ��Ǳ߽�ڵ�

***

## &nbsp; initial.h

**����һЩ��ʼ���ĺ��������������ʼ���ͳ���ֵ������**

### ����

+ ���ֱ�������  ***���Էŵ�config��ȥ***
  + �����������ڵ�����
  + �ٶȡ���������(centered)����Ҫ�������鱣������ʱ����ֵ��
  + һЩ��ʱ�õ��ı���������
    + �߳����ⷨ����������Щ��ô���棿����***���ñ��������ߵı߳���ֻ�ñ���ÿ���ڵ�������Ӧ�ı߳����ⷨ����ͬ��***��
    + �ڵ��ٶȣ����ܷ�ֱ�Ӷ��嵽�ṹ���У���
    + �ܶȣ�������ֻ����$t_n$ʱ�̵ĵ�ֵ���ó�ʼ������Jacobi����õ�����ֵ��������ɺ�������***�ܶ���ÿ��Gauss�ڵ������ڵ㶼��Ҫ������***��

### ����

+ ���񡢽ڵ��ʼ��
  + ��¼����ͽڵ��ţ�һά��źͶ�ά�ο����$(i,j)$��
  + ��¼�ڵ���������$(x,y)$
  + ��¼����Ķ��㣨����ĳ��˳��
  + ��¼�ڵ���������񣨰���ĳ��˳��
+ ��ֵ����
+ ��������
+ ����ÿ����Ԫ����������ÿ����Ԫ��������ͬ����Ҫÿ���������ڽṹ���У�
+ ��������ʼ��
  + ���񡢽ڵ��ʼ��
  + Jacobi
  + ���ġ���������
  + ������

***

## &nbsp; nodal_solver.h

**�ڵ���������������ڵ��ٶȺͽ���Ӧ����**

+ �ڵ��ٶ�
  + ��������ʶ��
  + ����
+ corner force

***

## &nbsp; time_evolution.h

**�����ʱ���ƽ��������ܶȡ��������ܡ���Ԫ�ٶȡ��ڵ��������ꡣ�����ܶ�ͨ��ǿ�����غ���㣬������ͨ����������RK���㡣**

+ �ܶȼ��㣬ǿ�����غ�
+ RK����ϵ���������($ \frac{dU^{(n+1)}}{dt} = L U^{(n)}$�е�$L$)
+ RK step1
+ RK step2

***

## &nbsp; limiter.h

**�Ա���������������������ʱ�䲽�������ơ�**

+ limiter
+ ʱ�䲽����ѡȡ

***

## &nbsp; plot.h

**�Խ����ͼ��**

***

## &nbsp; main

```flow
st=>start: Start
e=>end
initial=>operation: Initialization
end_time=>condition: t >=T?
EOS=>operation: EOS
limiter=>operation: limiter
RK=>operation: RK
ns=>operation: nodal solver

st->initial->EOS->limiter->ns->RK->end_time
end_time(yes)->e
end_time(no)->EOS
```