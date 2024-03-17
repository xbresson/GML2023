import torch
import time
import pickle
import dgl


# class of atom and bond dictionaries
class Dictionary:
    """
    word2idx and idx2word are mappings from words to idx and vice versa
    word2idx is a dictionary
    idx2word is a list
    word2num_occurence compute the number of times a given word has been added to the dictionary
    idx2num_occurence do the same, but with the index of the word rather than the word itself.
    """
    def __init__(self):
        self.word2idx = {}
        self.idx2word = []
        self.word2num_occurence = {}
        self.idx2num_occurence = []
    def add_word(self, word):
        if word not in self.word2idx:
            # dictionaries
            self.idx2word.append(word)
            self.word2idx[word] = len(self.idx2word) - 1
            # stats
            self.idx2num_occurence.append(0)
            self.word2num_occurence[word] = 0
        # increase counters    
        self.word2num_occurence[word]+=1
        self.idx2num_occurence[  self.word2idx[word]  ] += 1
    def get_rid_of_rare_words(self, min_num_occurence):
        new_idx2word = [ word for word in self.idx2word if self.word2num_occurence[word] >= min_num_occurence  ]
        new_word2idx = { word: idx  for idx,word in enumerate(new_idx2word) }         
        new_idx2num_occurence = [ self.word2num_occurence[word] for word in new_idx2word]   
        new_word2num_occurence = { word: self.word2num_occurence[word]  for word in new_idx2word } 
        self.word2idx = new_word2idx
        self.idx2word = new_idx2word
        self.word2num_occurence = new_word2num_occurence
        self.idx2num_occurence = new_idx2num_occurence
    def show(self):
        for idx, word in enumerate(self.idx2word):
            print(idx,'\t', word,'\t number of occurences = {}'.format(self.idx2num_occurence[idx]))
    def __len__(self):
        return len(self.idx2word)


class Molecule:
    """
    A molecule object contains the following attributes:
        ; molecule.num_atom : nb of atoms, an integer (N)
        ; molecule.atom_type : pytorch tensor of size N, each element is an atom type, an integer between 0 and num_atom_type-1
        ; molecule.atom_type_pe : pytorch tensor of size N, each element is an atom type positional encoding, an integer between 0 and num_atom-1
        ; molecule.bond_type : pytorch tensor of size N x N, each element is a bond type, an integer between 0 and num_bond_type-1 
        ; molecule.bag_of_atoms : pytorch tensor of size num_atom_type, histogram of atoms in the molecule
        ; molecule.logP_SA_cycle_normalized : the chemical property to regress, a pytorch float variable
        ; molecule.smile : the smile representation of the molecule for rdkit, a string   
    """
    def __init__(self, num_atom, num_atom_type):
        self.num_atom       = num_atom
        self.atom_type      = torch.zeros( num_atom , dtype=torch.long )
        self.atom_type_pe   = torch.zeros( num_atom , dtype=torch.long )
        self.bond_type      = torch.zeros( num_atom , num_atom, dtype=torch.long )
        self.bag_of_atoms   = torch.zeros( num_atom_type, dtype=torch.long)
        self.logP_SA        = torch.zeros( 1, dtype=torch.float)
        self.logP_SA_cycle_normalized  = torch.zeros( 1, dtype=torch.float)
        self.smile  = ''
    def set_bag_of_atoms(self):
        for tp in self.atom_type:
                self.bag_of_atoms[tp.item()] += 1
    def set_atom_type_pe(self):
        histogram={}
        for idx, tp in enumerate(self.atom_type):
            tpp=tp.item()
            if tpp not in histogram:
                histogram[tpp] = 0
            else:
                histogram[tpp] += 1
            self.atom_type_pe[idx] = histogram[tpp]
    def shuffle_indexing(self):
        idx = torch.randperm(self.num_atom)
        self.atom_type = self.atom_type[idx]
        self.atom_type_pe = self.atom_type_pe[idx]
        self.bond_type = self.bond_type[idx][:,idx]
        return idx
    def __len__(self):
        return self.num_atom

        
class MoleculeDGL(torch.utils.data.Dataset):
    def __init__(self, data_dir, split):
        self.split = split
        with open(data_dir + "/%s_pytorch.pkl" % split,"rb") as f:
            self.data = pickle.load(f)
        num_graphs = len(self.data)
        self.graph_lists = []
        self.graph_labels = []
        self.num_graphs = num_graphs
        self._prepare()
    def _prepare(self):
        print("preparing %d graphs for the %s set..." % (self.num_graphs, self.split.upper()))
        for molecule in self.data:
            node_features = molecule.atom_type.long()
            adj = molecule.bond_type
            edge_list = (adj != 0).nonzero() # converting adj matrix to edge_list
            edge_idxs_in_adj = edge_list.split(1, dim=1)
            edge_features = adj[edge_idxs_in_adj].reshape(-1).long()
            # version 2020 -- tmp
            #g = dgl.DGLGraph() # create the DGL graph
            #g.add_nodes(molecule.num_atom)
            #for src, dst in edge_list:
            #    g.add_edges(src.item(), dst.item())  
            # version 2023
            g = dgl.graph((edge_list[:,0], edge_list[:,1]), num_nodes=molecule.num_atom) # create the DGL graph  
            g.ndata['feat'] = node_features
            g.edata['feat'] = edge_features
            self.graph_lists.append(g)
            self.graph_labels.append(molecule.logP_SA_cycle_normalized)
    def __len__(self):
        return self.num_graphs
    def __getitem__(self, idx): # collate requires a method __getitem__ in the class
        return self.graph_lists[idx], self.graph_labels[idx]


class MoleculeDataset(torch.utils.data.Dataset):
    def __init__(self, data_name, data_dir):
        start = time.time()
        print("Loading datasets %s_dgl..." % (data_name))
        with open(data_dir + 'train_dgl.pkl',"rb") as f:
            self.train = pickle.load(f)
        with open(data_dir + 'val_dgl.pkl',"rb") as f:
            self.val = pickle.load(f)
        with open(data_dir + 'test_dgl.pkl',"rb") as f:
            self.test = pickle.load(f)
        print('train, test, val sizes :',len(self.train),len(self.test),len(self.val))
        print("Time: {:.4f}s".format(time.time()-start))
    # form a mini batch from a given list of samples = [(graph, label) pairs]
    # collate requires a method __getitem__ in the graph class used
    def collate(self, samples):
        # Input sample is a list of pairs (graph, label)
        graphs, labels = map(list, zip(*samples))
        batch_graphs = dgl.batch(graphs)
        batch_labels = torch.stack(labels)
        # Normalization w.r.t. graph sizes
        tab_sizes_n = [ graphs[i].number_of_nodes() for i in range(len(graphs))]
        tab_norm_n = [ torch.FloatTensor(size,1).fill_(1./float(size)) for size in tab_sizes_n ]
        batch_norm_n = torch.cat(tab_norm_n).sqrt()  
        tab_sizes_e = [ graphs[i].number_of_edges() for i in range(len(graphs))]
        tab_norm_e = [ torch.FloatTensor(size,1).fill_(1./float(size)) for size in tab_sizes_e ]
        batch_norm_e = torch.cat(tab_norm_e).sqrt()
        return batch_graphs, batch_labels, batch_norm_n, batch_norm_e
    

# Laplacian eigenvectors
def compute_LapEig(g, pos_enc_dim): # input g is a DGL graph
    Adj = g.adj().to_dense() # Adjacency matrix
    Dn = ( g.in_degrees()** -0.5 ).diag() # Inverse and sqrt of degree matrix
    Lap = torch.eye(g.number_of_nodes()) - Dn.matmul(Adj).matmul(Dn) # Laplacian operator
    EigVal, EigVec = torch.linalg.eig(Lap) # Compute full EVD
    EigVal, EigVec = EigVal.real, EigVec.real # make eig real
    EigVec = EigVec[:, EigVal.argsort()] # sort in increasing order of eigenvalues
    EigVec = EigVec[:,1:pos_enc_dim+1] # select the first non-trivial "pos_enc_dim" eigenvector
    return EigVec




