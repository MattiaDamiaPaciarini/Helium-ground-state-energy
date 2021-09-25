# Helium-ground-state-energy
The repository is composed by two different codes. Correct_Energy_Value.py uses the Metropolis algorithm and the VMC technique in order to find the correct value of the Helium nucleus ground state energy. Find_Optimal_Pars.py, instead, is able to find the optimal values of the variational parameters used in VMC technique starting in a neighbourhood of the optimal ones.
# Correct_Energy_Value.py Description
Starting from a random particles' configuration, Correct_Energy_Value.py uses the Metropolis algorithm to accept or not the following particles' positions, by comparing the wavefunction in the two different arrengements. In both cases, the algorithm computes the energy value depending on particles' positions. After a choosen number of iterations, the mean energy value is provided.
The technique used to compute the energy value is known as Variational MonteCarlo. It consists on evaluating the energy expectation value via a trial wavefunction <img src="https://render.githubusercontent.com/render/math?math=\psi(R,\alpha)"> depending on a set of parameters <img src="https://render.githubusercontent.com/render/math?math={\alpha_i}"> and suggested by system's symmetries. Then the parameters are tuned in order to reach the minimum of the energy. The code is already provided by the set of the optimal parameters.
Moreover, the code provides the user a collection of graphics, which allows a better understanding of particles' dynamics.
# Find_Optimal_Pars.py Description
Find_Optimal_Pars.py is a completion of the previous code. Starting from a set of parameters in a neighbourhood of the optimal one, it computes the average energy and finds the set of parameters which give the ground state energy. In order to optimize the algorithm, the so-called "Reweighting method" is employed: the energy value is evaluated by considering a different weight in front of it. To visualize the energy path followed to reach the minimum, a scatter plot is output
