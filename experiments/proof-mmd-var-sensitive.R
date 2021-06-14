devtools::load_all()

n = 1000


x1 = rnorm(n)
x2 = rnorm(n, sd = 1 + (x1 > 0) * 1)
df1 = data.frame(x1, x2)

x1a = rnorm(n)
x2a = rnorm(n, sd = 1 + (x1a > 0) * 1)
dfa = data.frame(x1 = x1a,x2 = x2a)


x1b = rnorm(n)
x2b = rnorm(n, sd = 1.5)
dfb = data.frame(x1 = x1b, x2 = x2b)

print(mmd2(df1, dfa))
print(mmd2(df1, dfb))


