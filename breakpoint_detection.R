# breakpoint detection
require(splines)

Month = c(-1.5,-0.5, 0.5 ,1.5, 2.5, 3.5, 4.5, 5.5, 6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5)
T = c(10.44444444, 5, 6.333333333, 7.333333333, 15.33333333,18.22222222,
      22.22222222,25.61111111, 29.72222222, 28.11111111,23.5,15.72222222, 12,6.111111111,
      5.222222222,6, 9.611111111)

model <- lm(T ~ ns(Month,7))
model
summary(model)

pdf("PNAS/plot/T_original.pdf", width=5, height=3.75)
X <- data.frame(Month = seq(-0.5,13, length=1000) ) # make an ordered sequence
Y <- predict(model,newdata = X) # calculate predictions for that sequence
plot(X$Month,Y,type="l", xlab = "Month", ylab = "Air temperature",
     xlim=c(0.1, 12.5), ylim=c(0, 30), xaxt="n") #check
axis(1, at = seq(1, 12, by = 1))
points(Month, T, col = "red", pch = 19)
dev.off()

pdf("PNAS/plot/T_Derivative.pdf", width=5, height=3.75)
dY <- diff(Y)/diff(X$Month)  # the derivative of your function
dX <- rowMeans(embed(X$Month,2)) # centers the X values for plotting
df = data.frame(dX, dY)
plot(dX,dY,type="l", xlab = "Month", ylab = "Month-to-month air temperature change",
     xlim=c(0, 12.5), ylim=c(-10,10), xaxt="n") #check
axis(1, at = seq(1, 12, by = 1))
points(c(2.374875, 9.271772), c(5.960660,  -7.325419), col = "red", pch = 19)
dev.off()
