# Load your data
cm_mouse <- get(load("data/cm_mouse.rda"))
cm_human <- get(load("data/cm_human.rda"))
cm_rat <- get(load("data/cm_rat.rda"))

CellCycleGenes_MOUSE <- get(load("data/CellCycleGenes_MOUSE.RData"))
CellCycleGenes_HUMAN <- get(load("data/CellCycleGenes_HUMAN.RData"))
CellCycleGenes_RAT <- get(load("data/CellCycleGenes_RAT.RData"))

# Combine cell cycle genes data
cellcyclegenes <- list(
  MOUSE = CellCycleGenes_MOUSE,
  HUMAN = CellCycleGenes_HUMAN,
  RAT = CellCycleGenes_RAT
)

# Save the data objects
usethis::use_data(cm_mouse, overwrite = TRUE)
usethis::use_data(cm_human, overwrite = TRUE)
usethis::use_data(cm_rat, overwrite = TRUE)
usethis::use_data(cellcyclegenes, overwrite = TRUE)
