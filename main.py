from parameters import *
from funk import *

def main():

    X, y = UploadDataset()
    X_train, X_test, y_train, y_test = OrganizeData(X, y)
    if optimize:
        MembraneOptimization(X_train, X_test, y_train, y_test)
    else:
        Membrane, genparam = ConstructMembraneFromDNA(upload_dna)
        Membrane, dna, success = AssessNetwork(Membrane, upload_dna, genparam, init=True)
        Membrane.SaveNetwork()
        print("SAVED MEMBRANE: " + str(Membrane.network_id))

if __name__ =="__main__":
    main()