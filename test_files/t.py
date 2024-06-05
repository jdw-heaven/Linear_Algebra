import torch

# 创建一个向量
x = torch.tensor([2, 3, 4, 5])

# 计算每个元素的倒数
y = torch.reciprocal(x)

# 打印结果
print("原始向量:", x)
print("倒数向量:", y)