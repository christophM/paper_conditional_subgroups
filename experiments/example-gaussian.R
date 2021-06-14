library("rpart")
library("MASS")
library("ggplot2")
library("patchwork")
library("dplyr")
devtools::load_all()

set.seed(42)
n = 100
sigm = matrix(c(1, 0.9, 0.9, 1), ncol = 2)


X = mvrnorm(n = n, mu = c(0, 0), Sigma = sigm)
X = data.frame(X)

rp = rpart(X2 ~ X1, data = X, control = rpart.control(maxdepth = 2))

splits = rp$splits[, "index"]
splits = data.frame(splits)
splits

p1 = ggplot(X) +
  scale_x_continuous(TeX("Feature x_1")) +
  geom_point(aes(x = X1, y = X2))

Xperm = X
Xperm$X2 = sample(Xperm$X2)
p2 = ggplot(Xperm) +
  geom_point(aes(x = X1, y = X2)) +
  scale_y_continuous("") +
  scale_x_continuous(TeX("Feature x_1"))
X$groups = factor(predict(rp, data = X))
Xperm2 = X %>% group_by(groups) %>% mutate(X2 = sample(X2))


p3 = ggplot(Xperm2) +
  geom_point(aes(x = X1, y = X2)) +
  geom_vline(aes(xintercept = splits), data = splits) +
  scale_x_continuous(TeX("Feature x_1")) +
  scale_y_continuous("")

pdf(file = sprintf("%s/gaussian.pdf", fig_dir), height = 3, width = 10)
p1 + p2 + p3
dev.off()
