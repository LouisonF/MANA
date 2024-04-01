# AIM : Assessing the robustness to the sampling approach

To reduce the computing time, we adapted the DEXOM approach.
The adapted DEXOM approach consist of :
* **Full Reaction-Enum procedure**
* **Stratified + Random sampling to select 1% of Reaction-Enum solutions**
* **Diversity-Enum starting from each selected solution**

# Assessment methodology

To assess how the solutions sampling may affect the results (*e.g.* the list of DARs), we will perform the adapted DEXOM approach several time (5) for the amiodarone usecase.
Practically it means:
    5 adapted DEXOM runs for 003016028014.CEL (amiodarone, ctrl, 24, sample1)
    5 adapted DEXOM runs for 003016028015.CEL (amiodarone, ctrl, 24, sample2)
    5 adapted DEXOM runs for 003016028020.CEL (amiodarone, high, 24, sample1)
    5 adapted DEXOM runs for 003016028021.CEL (amiodarone, high, 24, sample2)

Since in each adapted DEXOM runs, a random sampling is performed in each range defined by the stratified sampling step, we will assess how the random solution sampling might affect the results.

# Assessment results