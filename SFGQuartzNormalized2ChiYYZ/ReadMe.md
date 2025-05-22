# 用于将基于石英Normalized的SFG光谱，计算它的Chi-effective，以及计算SSP对应的YYZ

本程序基于Matlab开发。

程序中提供了两个例子，分别是水和苯的光谱。

## 使用方法：

- 将需要用到的介质的波长和折射率，新建csv表格存储在RefractiveIndex文件夹中，第一列设置为波长（单位um），第二列设置为折射率（n）
- 基于两个例子，对程序进行修改。
  - 修改实验构型，比如入射角，可见的波长与红外波长范围等
  - 修改路径依赖，如需要计算的光谱文件路径和计算折射率函数中使用的RefractiveIndex文件路径
  - 根据需要计算的内容修改程序


## 开发者信息
章力，西湖大学理学院化学系，王鸿飞实验室

zhangli@westlake.edu.cn


课题组主页https://hfwgroup.westlake.edu.cn/