library(gNanoPkg)

resultsRoot = "/Users/jcur002/Dropbox/Code/git/gNano/systematic/"

# g0.df = loadResults("g-0", resultsRoot = resultsRoot)
# g0.ll = calcLogLik(g0.df, "g")

ln0.df = loadResults("ln-0", resultsRoot = resultsRoot)
ln0.ll = calcLogLik(ln0.df, "ln")

ln1.df = loadResults("ln-1", resultsRoot = resultsRoot)
ln1.ll = calcLogLik(ln1.df, "ln")

ln2.df = loadResults("ln-2", resultsRoot = resultsRoot)
ln2.ll = calcLogLik(ln2.df, "ln")

ln3.df = loadResults("ln-3", resultsRoot = resultsRoot)
ln3.ll = calcLogLik(ln3.df, "ln")

ln4.df = loadResults("ln-4", resultsRoot = resultsRoot)
ln4.ll = calcLogLik(ln4.df, "ln")


compareLL(ln1.ll,ln2.ll,ln3.ll, ln4.ll, labels = c("ln-1", "ln-2", "ln-3", "ln-4"))


g1.df = loadResults("g-1", resultsRoot = resultsRoot)
g1.ll = calcLogLik(g1.df, "g")

compareLL(g0.ll, g1.ll, labesl = c("g-0", "g-1"))

sn1.df = loadResults("sn-1", resultsRoot = resultsRoot)
sn1.ll = calcLogLik(sn1.df, "sn")

sn2.df = loadResults("sn-2", resultsRoot = resultsRoot)
sn2.ll = calcLogLik(sn2.df, "sn")

sn3.df = loadResults("sn-3", resultsRoot = resultsRoot)
obsExpectPlot(sn3.df, "sn")
obsExpectPlot(sn3.df, "sn", TRUE)

sn3.ll = calcLogLik(sn3.df, "sn")

compareLL(sn2.ll, sn3.ll, labels = c("sn-2", "sn-3"))

sn4.df = loadResults("sn-4", resultsRoot = resultsRoot)
sn4.ll = calcLogLik(sn4.df, "sn")

compareLL(sn1.ll, sn2.ll, sn3.ll, sn4.ll, labels = c("sn-1", "sn-2", "sn-3", "sn-4"))

ln4.df = loadResults("ln-4", resultsRoot = resultsRoot)
ln4.ll = calcLogLik(ln4.df, "ln")
obsExpectPlot(ln4.df, "ln")
obsExpectPlot(ln4.df, "ln", TRUE)

compareLL(sn3.ll, ln4.ll, "sn-3", "ln-4")


sn3.df = loadResults("sn-3", resultsRoot = resultsRoot)
sn3.ll = calcLogLik(sn3.df, "sn")
sn3.plot = obsExpectPlot(sn3.df, "sn", TRUE)

ln4.df = loadResults("ln-4", resultsRoot = resultsRoot)
ln4.ll = calcLogLik(ln4.df, "ln")
ln4.plot = obsExpectPlot(ln4.df, "ln", TRUE)

compareLL(ln4.ll, sn3.ll, "ln-4", "sn-3")

p = sn3.plot$p +
  geom_point(data = ln4.plot$plot.df, aes(y = residuals, x = mean), col = "green")
p

g4.df = loadResults("g-4", resultsRoot = resultsRoot)
g4.plot = obsExpectPlot(g4.df, "g", TRUE)

p = g4.plot$p +
  geom_point(data = ln4.plot$plot.df, aes(y = residuals, x = mean), col = "green")
p
