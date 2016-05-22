# simulated random data
df <- data.frame(data = rnorm(60), expand.grid(time=1:10, group=1:2, id=1:3)) 
# apply test
res <- TPDT(df)
# print result
res
# plotting method
plot(res)
# illustration of concept
plot(res, plottype = 3)
