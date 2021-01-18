GECCO 2019 Competition - Bi-objective Traveling Thief Problem
============================================================

Please have a look at the competition outline `Here 
<https://www.egr.msu.edu/coinlab/blankjul/gecco19-thief/>`_.

Requirements
------------------------------------------------------------
- Java 8
- (Maven)

Installation
------------------------------------------------------------

Feel free to use the IDE of your choice and import the Maven Project.


Structure
------------------------------------------------------------

In the following the project structure is explained:

::

    gecco19-thief
    ├── Runner.java: Execute an algorithm on all competition instance and to save the file in the derired format.
    ├── Competition.java: Contains the instance names to be solved and the maximum limit of solutions to submit.
    ├── model
        ├── TravelingThiefProblem.java: The problem object used to evaluate the problem for a given tour and packing plan.
        ├── Solution.java: Object to store the results of the evaluate function.
        └── NonDominatedSet.java: Example implementation of a non-dominated set. Can be done faster/better.
    ├── algorithms
        ├── Algorithm: Interface for the algorithm to be implemented from.
        ├── NtgaAlgorithm: NTGA algorithm implemented by Team Lambda.
        └── LambdaAlgorithm: NSGA II algorithm implemented by Team Lambda.



Getting Started
------------------------------------------------------------

Please have a look at Team Lambda's implementation of the NSGA II algorithm in the file LambdaAlgorithm.java. The results of all the problem sets can be found in the folder /gecco19-thief/results.
To check the algorithm against a particular problem please change the instanceToRun variable in the Runner.java file and run. 
More experimentation can done by changing the tournamentSize, mutationFactor, clonePrevention and totalGenerations variables in the LambdaAlgorithm.


