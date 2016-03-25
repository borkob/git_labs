ff <- read.table("walkLength.txt")
n <- length(ff$V1)
s <- sd(ff$V1)
a <- mean(ff$V1)
error <- qnorm(0.975)*s/sqrt(n)
left <- a-error
right <- a+error
print(paste("mean walk length:",a))
print(paste("st of walk length:",s))
print(paste("runs:",n))
print(paste("walk length confidence interval:",left,right))
