2016-07-08
1.全局LLE效果极差
2.IID_adapted为我重写作者Junho的代码
3.多数函数放在utils中
4.使用demo进行操作

2016-07-11
Add IID_slic.m 
* IID_slic(I, I, depth, 1) means not use filter 
Add demo-NYU.m
Add demo-Sintel.m
* add option such as select IID_slic or not, select cropping(for NYU dataset) or not
Improve Local Retinex Constraints Computation
* LocalNormalConstraintMatrix.m
* LocalReflectanceConstraintMatrix.m