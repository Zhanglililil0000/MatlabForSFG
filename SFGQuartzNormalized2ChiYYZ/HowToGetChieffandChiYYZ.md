# 本程序基于Matlab开发

# 液气界面 HR-BB-SFG-VS 结果归一化与结果分析

​	飞秒宽带SFG-VS光谱，是通过窄带皮秒可见光和飞秒宽带红外光共同作用的过程。由于飞秒宽带红外光在不同波长的能量分布存在差异，所以得到的SFG-VS结果需要通过归一化的方式，获得准确的光谱形状。常见方法是使用强度稳定的介质，比如金、GaAs晶体、石英晶体作为参考，对得到的光谱结果进行归一化。在实验中，SFG-VS的强度可以写作
$$
I=\frac{8\pi ^3\omega ^2\sec ^2\beta}{c^3n_I\left( \omega \right) n_I\left( \omega _1 \right) n_I\left( \omega _2 \right)}\left| \chi _{eff}^{\left( 2 \right)} \right|^2I\left( \omega _1 \right) I\left( \omega _2 \right)
\tag{1}
$$
​	那么归一化过程就可以写作
$$
I_{Normalized}=\frac{I_s}{I_Q}=\frac{\left| \chi _{eff,s}^{\left( 2 \right)} \right|^2}{\left| \chi _{eff,Q}^{\left( 2 \right)} \right|^2}
\tag{2}
$$
​	其中 $I_{s}$ 是样品的SFG信号， $I_{Q}$ 是相同条件下石英的SFG信号。式(1)中的光速、入射相折射率、SFG折射角度等这些各项物理参数，以及入射光能量 $I(\omega_1), I(\omega_2)$ 在归一化过程中被消除掉了，这个时候得到的就是对石英归一化的SFG结果 $I_{Normalized}$ 。所以一旦已知该条件下的石英的 $\left| \chi _{eff,Q}^{\left( 2 \right)} \right|^2$ ，就可以计算出样品的 $\left| \chi _{eff,s}^{\left( 2 \right)} \right|^2$
$$
\left| \chi _{eff,s}^{\left( 2 \right)} \right|^2=I_{Normalized}\times \left| \chi _{eff,Q}^{\left( 2 \right)} \right|^2
\tag{3}
$$
​	在SSP偏振组合下，$C_{\infin,v}$ 的宏观对称性条件下，$\chi_{eff}$ 就主要由一项二阶非线性极化张量的元贡献
$$
\chi_{ssp}=L_{yy}(\omega)L_{yy}(\omega_{1})L_{zz}(\omega_{2})sin\beta \cdot \chi_{yyz}
$$
​	其中L是菲涅尔因子，所以为了计算 $\chi_{yyz}$ ，需要计算得到菲涅尔因子L，然后把结果修正回去。

​	那么这里的修正方法就很直观了，先得到 Normalized 的 SFG-VS 结果，然后分别乘 $\left| \chi _{eff,Q}^{\left( 2 \right)} \right|^2$ 和除去菲涅尔因子和入射角项，就可以得到最终的 $\chi_{yyz}$ ，也是我们最终在拟合光谱时候应该拟合的结果。



## 附录A 石英的折射率计算

