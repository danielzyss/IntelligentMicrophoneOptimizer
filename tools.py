from parameters import *
from Membrane import MembraneSystem
from Particle import Particle
from Spring import Spring

def CreateNetworkElements(genparam):

    if net_dim==2:
        x, y = np.linspace(0, genparam['xlength'], genparam['row']), np.linspace(0, genparam['ylength'], genparam['col'])
        xx, yy = np.meshgrid(x, y)
        grid_points = np.stack((xx.flatten(), yy.flatten()), axis=-1)
        sigma_pos = ((genparam['xlength'] + genparam['ylength']) / 2) / (2 * (genparam['col'] +genparam['row']) / 2)
    else:
        raise NameError('Incorrect Network Initialization Dimensions (net_dim must be 2)')

    feed_gen = np.random.uniform(2, grid_points.shape[0])
    feed_id = np.random.choice(np.arange(0, grid_points.shape[0]), round(feed_gen), replace=False).tolist()

    Particles = []
    for i_p, p in enumerate(grid_points):
        fixed = IsFixed(p, genparam)
        input = IsInput(p, fixed, genparam)

        if not fixed and not input:
            p= np.abs(p+ np.random.normal(0, sigma_pos, net_dim))

        if i_p in feed_id:
            if not fixed and not input:
                # w_feed = np.random.normal(genparam['feedgen'], genparam['sig']) * genparam['overfeed']
                w_feed = np.random.normal(0.0, 1.0) * genparam['overfeed']
            else:
                feed_id[feed_id.index(i_p)] = i_p+1
                w_feed = 0.0
        else:
            w_feed = 0.0

        if input:
            # w_input = np.random.normal(genparam['ingen'], genparam['sig']) * genparam['overin']
            w_input = np.random.normal(0.0, 1.0) *genparam['overin']
        else:
            w_input = 0.0

        # m = np.random.normal(genparam['massgen'], genparam['sig'])
        m = 1.0
        Particles.append(Particle(i_p, p, m, fixed, w_feed, w_input))
        grid_points[i_p] = p
    Particles = np.array(Particles)

    if genparam['del']:
        mesh = Delaunay(grid_points)
    else:
        mesh = RandomConnectionMesh(grid_points, genparam)

    edges, nb_edges = BuildEdges(mesh, genparam['del'])

    Springs = []

    for j, e_j in enumerate(edges):
        l0 = euclidean(grid_points[e_j[0]], grid_points[e_j[1]])

        if quadratic_spring:
            k = np.abs(np.random.normal(genparam['stiffgen'], genparam['sig'], 2))
            d = np.abs(np.random.normal(genparam['dampgen'], genparam['sig'], 2))
        else:
            k = np.abs(np.random.normal(genparam['stiffgen'], genparam['sig']))
            d = np.abs(np.random.normal(genparam['dampgen'], genparam['sig']))

        new_spring = Spring(l0, e_j[0], e_j[1], k, d)
        Springs.append(new_spring)
    Springs = np.array(Springs)

    return Particles, Springs, mesh, nb_edges

def BuildEdges(tri, delo):

    edges = []
    visitedEdges = []

    if delo:
        for s in tri.simplices:
            for e in list(itertools.combinations(s, 2)):
                if e not in visitedEdges:
                    edges.append(e)
                    visitedEdges.append((e[0], e[1]))
                    visitedEdges.append((e[1], e[0]))
    else:
        for i in range(0, tri.shape[0]):
            for j in tri[i]:
                if i != j and [i, j] not in visitedEdges and [j, i] not in visitedEdges:
                    edges.append(np.array([i, j]))
                    visitedEdges.append([i, j])

    edges = np.array(edges)
    nbOfEdges = edges.shape[0]

    return edges, nbOfEdges

def IsFixed(pt, gennparam):

    if net_dim==2:
        if pt[1] == gennparam['ylength']:
            return True
        elif pt[0] == gennparam['xlength']:
            return True
        elif pt[1] == 0:
            return True
        else:
            return False

def IsInput(pt, fixed, genparam):

    if net_dim==2:
        if pt[0] == 0 and pt[1] != 0 and pt[1] != genparam['ylength']:
            return True
        else:
            return False

def RandomConnectionMesh(points, genparam):

    Neigh = NearestNeighbors()
    Neigh.fit(points)
    mesh = []
    for pt in points:
        r_norm = np.random.normal(0, ((genparam['xlength']+genparam['ylength'])/2)/((genparam['row']+genparam['col'])/2))
        dist, neighbours_idx = Neigh.radius_neighbors(pt.reshape((1, net_dim)), np.abs(r_norm))
        mesh.append(neighbours_idx[0])

    return np.array(mesh)

def AssessNetwork(Membrane, dna, genparam, init=False):
    n_it =0
    step = np.hstack((np.zeros(20), np.ones(washout_max_permitted_time/dt - 20)))/genparam['overin']
    M = np.zeros((Membrane.nb_edge, int(washout_max_permitted_time/dt)))

    validated = False

    while not validated:
        if n_it>max_children_gen_it:
            return Membrane, dna, False

        no_nan = False
        while not no_nan:
            if n_it > max_children_gen_it:
                return Membrane, dna, False

            if cplusplus:
                M = RunWithCplusplus(Membrane, 'loopless', M=M, force=step)
                if not np.isnan(M).any():
                    no_nan = True
                else:
                    if init:
                        dna = "".join(np.random.choice(["0", "1"], size=(len_dna)).tolist())
                        Membrane, genparam = ConstructMembraneFromDNA(dna)
                        M = np.zeros((Membrane.nb_edge, int(washout_max_permitted_time / dt)))
                        step = np.hstack((np.zeros(20), np.ones(washout_max_permitted_time / dt - 20))) / genparam[
                            'overin']
                    else:
                        n_it+=1
                        Membrane, genparam= ConstructMembraneFromDNA(dna)
                        M = np.zeros((Membrane.nb_edge, int(washout_max_permitted_time / dt)))
                        step = np.hstack((np.zeros(20), np.ones(washout_max_permitted_time / dt - 20))) / genparam[
                            'overin']

            else:
                try:
                    M = Membrane.RunLoopless(M, step)
                    no_nan = True
                except:
                    if init:
                        dna = "".join(np.random.choice(["0", "1"], size=(len_dna)).tolist())
                        Membrane, genparam= ConstructMembraneFromDNA(dna)
                        M = np.zeros((Membrane.nb_edge, int(washout_max_permitted_time / dt)))
                        step = np.hstack((np.zeros(20), np.ones(washout_max_permitted_time / dt - 20))) / genparam[
                            'overin']
                    else:
                        n_it += 1
                        Membrane, genparam = ConstructMembraneFromDNA(dna)
                        M = np.zeros((Membrane.nb_edge, int(washout_max_permitted_time / dt)))
                        step = np.hstack((np.zeros(20), np.ones(washout_max_permitted_time / dt - 20))) / genparam[
                            'overin']

        grad_m = []
        for m in M:
            delta = np.abs(np.gradient(m, dt))
            grad_m.append(delta/np.sqrt(genparam['xlength']**2+ genparam['ylength']**2))

        std_m = np.std(grad_m, axis=0)
        std_m_norm = (std_m-std_m.min())/(std_m.max()-std_m.min())
        washout_idx = np.where(std_m_norm > washout_criteria)[0]

        if washout_idx.shape[0]>0 and np.max(grad_m)<1:
            washout_idx = max(washout_idx)
            if washout_idx<(washout_max_permitted_time/dt - 1):
                Membrane.washout_time = dt*washout_idx
                Membrane.InitialState = Membrane.ConstructStruct(Membrane.Particles, Membrane.Springs,
                                                                 Membrane.Mesh, Membrane.nb_edge,
                                                                 Membrane.network_id)
                return Membrane, dna, True
            else:
                if init:
                    dna = "".join(np.random.choice(["0", "1"], size=(len_dna)).tolist())
                    Membrane, genparam= ConstructMembraneFromDNA(dna)
                    M = np.zeros((Membrane.nb_edge, int(washout_max_permitted_time / dt)))
                    step = np.hstack((np.zeros(20), np.ones(washout_max_permitted_time / dt - 20))) / genparam['overin']
                else:
                    n_it += 1
                    Membrane, genparam= ConstructMembraneFromDNA(dna)
                    M = np.zeros((Membrane.nb_edge, int(washout_max_permitted_time / dt)))
                    step = np.hstack((np.zeros(20), np.ones(washout_max_permitted_time / dt - 20))) / genparam['overin']
        else:
            if init:
                dna = "".join(np.random.choice(["0", "1"], size=(len_dna)).tolist())
                Membrane, genparam= ConstructMembraneFromDNA(dna)
                M = np.zeros((Membrane.nb_edge, int(washout_max_permitted_time / dt)))
                step = np.hstack((np.zeros(20), np.ones(washout_max_permitted_time / dt - 20))) / genparam['overin']
            else:
                n_it += 1
                Membrane, genparam= ConstructMembraneFromDNA(dna)
                M = np.zeros((Membrane.nb_edge, int(washout_max_permitted_time / dt)))
                step = np.hstack((np.zeros(20), np.ones(washout_max_permitted_time / dt - 20))) / genparam['overin']

def LoadMembrane(fileid=membrane_id):
    InitialState = pickle.load(open('Membranes/' + fileid, 'rb'))

    particles = []
    springs = []

    for i_p in range(0, len(InitialState['particles']['pos'])):
        p = Particle(i_p,
                     InitialState['particles']['pos'][i_p],
                     InitialState['particles']['mass'][i_p],
                     InitialState['particles']['fixed'][i_p],
                     InitialState['particles']['w_feed'][i_p],
                     InitialState['particles']['w_in'][i_p])
        particles.append(p)

    for i_s in range(0, len(InitialState['springs']['l0'])):
        s = Spring(InitialState['springs']['l0'][i_s],
                   InitialState['springs']['p1'][i_s],
                   InitialState['springs']['p2'][i_s],
                   InitialState['springs']['k'][i_s],
                   InitialState['springs']['d'][i_s])
        springs.append(s)

    m = InitialState['mesh']
    n = InitialState['nb_edge']
    id = InitialState['id']
    w = InitialState['washout_time']
    dim = InitialState['net_dim']
    dna = InitialState['dna']

    M = MembraneSystem(particles, springs, m, n, dim, dna)
    M.id = id
    M.washout_time = w

    return M

def OrganizeData(X, y):
    X = X[:nb_classes]
    learning_feed = []
    for i, c in enumerate(X):
        learning_feed.append([])
        for x in c:
            learning_feed[i].append(np.full(x.shape[0], i, dtype=np.int))
            # (i + 1) * ((-1) ** (i + 1))
    learning_feed = np.array(learning_feed)

    Yset = learning_feed.reshape((learning_feed.shape[0] * learning_feed.shape[1]))
    Xset = X.reshape((X.shape[0] * X.shape[1]))

    X_train, X_test, y_train, y_test = train_test_split(Xset, Yset, train_size=0.8)

    return X_train, X_test, y_train, y_test

def UploadDataset():
    X = np.load('Data/voweldataset/sounds.npy')
    y = np.load('Data/voweldataset/target.npy')
    return X,y

def PlotTraining(Membrane, M_train, w, Y_train, vowel_length):

    regressed = np.matmul(w, M_train)

    regressed_r = np.array([])
    for i in range(0, train_batch_size):
        mean = np.round(np.mean(regressed[i*vowel_length:i*vowel_length+vowel_length]))
        regressed_r = np.append(regressed_r, np.full((vowel_length), mean))

    plt.plot(np.matmul(w, M_train), 'b', label='Linear Regression')
    plt.plot(Y_train, 'r', label='Expected Value')
    plt.plot(regressed_r, '--g', label='Running Average')
    plt.ylabel('Class')
    plt.xlabel('Time')
    plt.title('Training Regression vs. Expected Classes')
    plt.legend()
    ftitle = 'Training-' + str(nb_classes) + 'classes-' + str(train_batch_size) + 'size-' + str(Membrane.network_id) + '.png'
    plt.savefig('figures/' + ftitle)
    plt.close()

    return regressed, regressed_r

def PlotTesting(Membrane, M_test, w, Y_test, vowel_length):

    predicted = np.matmul(w, M_test)
    predicted_r = np.array([])
    for i in range(0, test_batch_size):
        mean = np.round(np.mean(predicted[i * vowel_length:i * vowel_length + vowel_length]))
        predicted_r = np.append(predicted_r, np.full((vowel_length), mean))

    plt.plot(predicted, 'b', label='Linear Regression')
    plt.plot(Y_test, 'r', label='Expected Value')
    plt.plot(predicted_r, '--g', label='Rounded Regression', )
    plt.ylabel('Class')
    plt.xlabel('Time')
    plt.title('Testing Regression vs. Expected Classes')
    plt.legend()
    ftitle = 'Testing-' + str(nb_classes) + 'classes-' + str(test_batch_size) + 'size-' + str(Membrane.network_id) + '.png'
    plt.savefig('figures/' + ftitle)
    plt.close()

    return predicted, predicted_r

def RunWithCplusplus(Membrane, action='openloop', M=np.array([]), force=np.array([]), y=np.array([]), w=np.array([])):

    command = "./fastmembrane "
    command += action + " "
    command += str(len(Membrane.Particles)) + " "
    command += str(len(Membrane.Springs)) + " "
    command += str(Membrane.washout_time) + " "
    command += "euler" + " "
    command += str(M.shape[0]) + " "
    command += str(M.shape[1]) + " "
    command += str(dt) + " "
    command += str(M.shape[1]*dt) + " "
    command += str(Membrane.net_dim) + " "

    if action=='openloop':
        np.savetxt("Cplusplus/openloop_force.txt", force, delimiter="\n", fmt='%1.6f')
        np.savetxt("Cplusplus/openloop_y.txt", y, delimiter="\n", fmt='%1.6f')
        command+= "openloop_force.txt" + " "
        command+= "openloop_y.txt" + " "
    elif action=='loopless':
        np.savetxt("Cplusplus/loopless_force.txt", force, delimiter="\n", fmt='%1.6f')
        command += "loopless_force.txt" + " "
    elif action=='closedloop':
        np.savetxt("Cplusplus/closedloop_force.txt", force, delimiter="\n", fmt='%1.6f')
        np.savetxt("Cplusplus/closedloop_w.txt", w, delimiter="\n", fmt='%1.6f')
        command += "closedloop_force.txt" + " "
        command += "closedloop_w.txt" + " "

    for p in Membrane.Particles:
        command+= str(p.pos[0]) + " "
    for p in Membrane.Particles:
        command+= str(p.pos[1]) + " "
    for p in Membrane.Particles:
        command+= str(p.mass) + " "
    for p in Membrane.Particles:
        command+= str(p.fixed) + " "
    for p in Membrane.Particles:
        command+= str(p.w_feed) + " "
    for p in Membrane.Particles:
        command+= str(p.w_in) + " "
    for s in Membrane.Springs:
        command+= str(s.l0) + " "
    for s in Membrane.Springs:
        command+= str(s.p1) + " "
    for s in Membrane.Springs:
        command+= str(s.p2) + " "
    for s in Membrane.Springs:
        command+= str(s.k1) + " "
    for s in Membrane.Springs:
        command+= str(s.k2) + " "
    for s in Membrane.Springs:
        command+= str(s.d1) + " "
    for s in Membrane.Springs:
        command+= str(s.d2) + " "
    if Membrane.net_dim==3:
        for p in Membrane.Particles:
            command+= str(p.pos[2]) + " "


    owd = os.getcwd()
    os.chdir('Cplusplus')


    try:
        out = subprocess.check_output(command.split())
        out = out[:len(out) - 1]
        M = np.array([float(i) for i in out.decode("utf-8").split("*")]).reshape(M.shape)
    except:
        os.chdir(owd)
        if action == 'openloop':
            os.remove("Cplusplus/openloop_force.txt")
            os.remove("Cplusplus/openloop_y.txt")
        elif action == 'loopless':
            os.remove("Cplusplus/loopless_force.txt")
        elif action == 'closedloop':
            os.remove("Cplusplus/closedloop_force.txt")
            os.remove("Cplusplus/closedloop_w.txt")
        return np.full(M.shape, np.nan)

    os.chdir(owd)

    if action=='openloop':
        os.remove("Cplusplus/openloop_force.txt")
        os.remove("Cplusplus/openloop_y.txt")
    elif action=='loopless':
        os.remove("Cplusplus/loopless_force.txt")
    elif action=='closedloop':
        os.remove("Cplusplus/closedloop_force.txt")
        os.remove("Cplusplus/closedloop_w.txt")


    return M

def PlotInitialNetwork(grid_points, edges, Particles ):
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for p in grid_points:
        f = IsFixed(p)
        if f:
            ax.scatter(p[0], p[1], p[2], c='r')
        elif IsInput(p, f):
            ax.scatter(p[0], p[1], p[2], c='g')
        else:
            ax.scatter(p[0], p[1], p[2], c='b')

    for e in edges:
        ax.plot([Particles[e[0]].pos[0], Particles[e[1]].pos[0]], [Particles[e[0]].pos[1], Particles[e[1]].pos[1]],
                [Particles[e[0]].pos[2], Particles[e[1]].pos[2]], 'b')

    plt.show()

#EVOLUTIONNARY:

def ConstructMembraneFromDNA(dna):

    generation_param ={}
    for i, gene in enumerate(gene_names):
        if i<2:
            generation_param[gene] = min(maxminval[i][0],
                                         int(dna[gene_cuts[i][0]:gene_cuts[i][1]], 2)+ maxminval[i][1])
        elif i>=2 and i<4:
            generation_param[gene] = max(maxminval[i],
                                         10**int(dna[gene_cuts[i][0]:gene_cuts[i][1]], 2))
        elif i==4:
            generation_param[gene] = bool(int(dna[gene_cuts[i]]))
        elif i>4:
            generation_param[gene] = 10**([-1,1][int(dna[gene_cuts[i][0]])]*int(dna[gene_cuts[i][0]+1:gene_cuts[i][1]],2))

    Particles, Springs, mesh, nb_edges = CreateNetworkElements(generation_param)
    Membrane = MembraneSystem(Particles, Springs, mesh, nb_edges, net_dim, dna)

    return Membrane, generation_param

def GenerateInitialPopulation(n):
    pop_id = []
    pop_dna = []
    pop_M = []
    for i in tqdm.tqdm(range(0, n), desc='Initial Population'):
        dna = "".join(np.random.choice(["0", "1"], size=(len_dna)).tolist())
        Membrane, genparam = ConstructMembraneFromDNA(dna)
        Membrane, dna, success = AssessNetwork(Membrane, dna, genparam, init=True)
        pop_dna.append(dna)

    return np.array(pop_dna)

def FitnessMembrane(dna, X_train, X_test, y_train, y_test):

    testing_accuracy = 0
    n=0
    idx_train = np.arange(X_train.shape[0])
    idx_test = np.arange(X_test.shape[0])

    while n<score_attempt:
        if n>2*score_attempt:
            return 999999.9999
        Membrane, genparam = ConstructMembraneFromDNA(dna)
        Membrane, dna, success = AssessNetwork(Membrane, dna, genparam)
        Membrane.InitialState = Membrane.ConstructStruct(Membrane.Particles,Membrane.Springs, Membrane.Mesh,Membrane.nb_edge,Membrane.network_id)
        np.random.shuffle(idx_train)
        np.random.shuffle(idx_test)

        X_train = X_train[idx_train]
        y_train = y_train[idx_train]
        X_test = X_test[idx_test]
        y_test = y_test[idx_test]

        M_train = []
        Y_train = []

        for b in range(0, train_batch_size):
            Membrane.resetNetwork()
            xb = X_train[b]
            yb = y_train[b]
            Mb = np.zeros((Membrane.nb_edge, yb.shape[0]))

            if cplusplus:
                Mb = RunWithCplusplus(Membrane, 'openloop', M=Mb, force=xb, y=yb)
            else:
                Mb = Membrane.RunOpenLoop(Mb, xb, yb)

            M_train.append(Mb[:, int(Membrane.washout_time / dt):])
            Y_train.append(y_train[b][int(Membrane.washout_time /dt):])

        M_train = np.hstack(M_train)
        M_train += np.random.normal(0.0, 0.001, M_train.shape)
        Y_train = np.hstack(Y_train)

        if np.isnan(M_train).any():
            continue

        M_train = sm.add_constant(M_train.T, has_constant='add').T

        model = sm.OLS(Y_train, M_train.T)
        results = model.fit()
        w = results.params

        vowel_length = []
        M_test = []
        Y_test = []
        for b in range(0, test_batch_size):
            Membrane.resetNetwork()
            xb = X_test[b]
            Mb = np.zeros((Membrane.nb_edge, xb.shape[0]))
            if cplusplus:
                Mb = RunWithCplusplus(Membrane, 'closedloop', M=Mb, force=xb, w=w)
            else:
                Mb = Membrane.RunClosedLoop(Mb, xb, w)
            M_test.append(Mb[:, int(Membrane.washout_time /dt):])
            Y_test.append(y_test[b][int(Membrane.washout_time /dt):])
            vowel_length.append(xb.shape[0])

        M_test = np.hstack(M_test)
        Y_test = np.hstack(Y_test)
        vowel_length = np.array(vowel_length)

        M_test = sm.add_constant(M_test.T, has_constant='add').T

        predicted = np.matmul(w, M_test)

        testing_accuracy += ComputeLogLoss(Y_test, predicted, vowel_length)
        n+=1


    return testing_accuracy/score_attempt

def ComputeLogLoss(Y_test, predicted, vowel_length):
    try:
        predicted_prob = []
        predicted_truth = []
        space = 0
        for b in range(0, test_batch_size):
            sample = predicted[space:space + vowel_length[b]]
            sample_p1 = round(np.mean(1 - sample))
            sample_p2 = round(np.mean(sample))
            predicted_prob.append(np.array([sample_p1, sample_p2]))
            predicted_truth.append(Y_test[space + 2])
            space += vowel_length[b]
        predicted_prob = np.array(predicted_prob)
        predicted_truth = np.array(predicted_truth)
        return sklearn.metrics.log_loss(predicted_truth, predicted_prob, labels=[0, 1])
    except:
        return 99999.9





def GenerateNewPopulation(Pop_dna, rank):

    #Conservation of elit
    NewPop_dna=Pop_dna[rank[:n_elit]].tolist()

    #Reproduction
    while len(NewPop_dna)<popsize:
        [p_a,p_b] = np.random.choice(rank[:crossover_size],2, replace=False)
        children_a, children_b = CrossOver(Pop_dna[p_a], Pop_dna[p_b])
        children_a = Mutation(children_a)
        children_b = Mutation(children_b)

        Memb_a, genparam_a = ConstructMembraneFromDNA(children_a)
        Memb_b, genparam_b = ConstructMembraneFromDNA(children_b)
        Memb_a, children_a, success_a = AssessNetwork(Memb_a,children_a, genparam_a)
        Memb_b, children_b, success_b = AssessNetwork(Memb_b, children_b, genparam_b)

        if success_a:
            NewPop_dna.append(children_a)
        if success_b:
            NewPop_dna.append(children_b)

    NewPop_dna = np.array(NewPop_dna)

    return NewPop_dna

def CrossOver(parent_a, parent_b):
    parent_a = list(parent_a)
    parent_b = list(parent_b)

    n_cuts = np.random.randint(0, np.ceil(len(gene_cuts)/2 - 1))*2+1
    cut_idx = np.random.choice(gene_cuts, n_cuts, replace=False)
    children_a = parent_a
    children_b = parent_b
    for jdx in cut_idx:
        if type(jdx)==int:
            children_a[jdx]= parent_b[jdx]
            children_b[jdx]= parent_a[jdx]
        else:
            children_a[jdx[0]:jdx[1]] = parent_b[jdx[0]:jdx[1]]
            children_b[jdx[0]:jdx[1]] = parent_a[jdx[0]:jdx[1]]

    children_a = "".join(children_a)
    children_b = "".join(children_b)

    return children_a, children_b

def Mutation(dna):
    dna = list(dna)
    #mutation 1:
    if np.random.uniform(0,1)<p_mut1:
        swapidx = np.random.randint(0, len(dna)-l_seg, n_seg)
        i= 0
        for idx in swapidx[::2]:
            a = dna[idx:idx+l_seg]
            b = dna[swapidx[i+1]:swapidx[i+1]+l_seg]
            dna[idx:idx+l_seg] = b
            dna[swapidx[i+1]:swapidx[i+1]+l_seg] = a
            i+=2

    #mutation 2:
    for i, b in enumerate(dna):
        if np.random.uniform(0,1)<p_mut2:
            if b=='0':
                dna[i]='1'
            else:
                dna[i]='0'

    return "".join(dna)
