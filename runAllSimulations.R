library(gNanoPkg)

debug(runSim)

data.df = readData()
simPath = here("systematic")

runSim(y~1, data = data.df, simPath = simPath, simRoot = "g-0")
runSim(y~profile, data = data.df, simPath = simPath, simRoot = "g-1")
runSim(y~profile + locus, data = data.df, simPath = simPath, simRoot = "g-2")
runSim(y~profile + locus + dye + X, data = data.df, simPath = simPath, simRoot = "g-3", genInits = TRUE)
runSim(y~profile + locus + dye + X + V, data = data.df, simPath = simPath, simRoot = "g-4",)

runSim(y~1, data = data.df, simPath = simPath, simRoot = "ln-0", responseDist = "normal")
runSim(y~profile, data = data.df, simPath = simPath, simRoot = "ln-1", responseDist = "normal")
runSim(y~profile + locus, data = data.df, simPath = simPath, simRoot = "ln-2", responseDist = "normal")
runSim(y~profile + locus + dye + X, data = data.df, simPath = simPath, simRoot = "ln-3", responseDist = "normal")
runSim(y~profile + locus + dye + X + V, data = data.df, simPath = simPath, simRoot = "ln-4", responseDist = "normal")

