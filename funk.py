from tools import *
from parameters import *

def GenerateMembrane():

    Particles, Springs, mesh, nb_edges = CreateNetworkElements()
    Membrane = AssessNetwork(MembraneSystem(Particles,Springs, mesh, nb_edges, net_dim))
    Membrane.SaveNetwork()

    return Membrane
def MembraneOptimization(X_train, X_test, y_train, y_test):

    Pop_dna = GenerateInitialPopulation(popsize)

    Pop_score = []
    for d in tqdm.tqdm(Pop_dna, desc="Fitness Generation 0"):
        score = FitnessMembrane(d, X_train, X_test, y_train, y_test)
        Pop_score.append(score)


    Pop_score = np.array(Pop_score)
    rank = np.argsort(Pop_score)

    for g in range(1, n_gen):
        
        print("Generation #" + str(g)+ " max fitness : " + str(Pop_score[rank[0]]))
        print("Elit Scores: ", Pop_score[rank[:n_elit]])
        for d in Pop_dna[rank[:n_elit]]: 
            print('Elit DNA: ' + d)
        
        Pop_dna = GenerateNewPopulation(Pop_dna, rank)

        Pop_score = []

        for d in tqdm.tqdm(Pop_dna, desc="Fitness Calculation " + str(g)):
            score = FitnessMembrane(d, X_train, X_test, y_train, y_test)
            Pop_score.append(score)

        Pop_score = np.array(Pop_score)
        rank = np.argsort(Pop_score)

    for i, d in enumerate(Pop_dna[rank[:10]]):
        print('Elit DNA #' + str(i)+ " " + d)
