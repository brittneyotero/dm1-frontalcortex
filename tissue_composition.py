import sys, os
import numpy as np


#Read in expression table of scRNA-seq data
#Combine all of the cells of the same type and average them
#Make a new expression table with average cell type expression
def cellTypeAvg(exprTable, outname):

        idx = {}
        expr = {}
        for line in open(exprTable):
                if line.startswith('Gene'):
                        header = line.strip().split('\t')
                        for i in range(len(header)):
                                if header[i] != 'Gene':
                                        cell = header[i].strip().split('_')[1]
                                        if cell != 'fetal' and cell != 'hybrid':
                                                if cell not in idx:
                                                        idx[cell] = [i]
                                                else:
                                                        idx[cell].append(i)
                else:
                        vals = line.strip().split('\t')
                        gene = vals[0]
                        expr[gene] = {}
                        for ctype in idx:
                                if ctype not in expr[gene]:
                                        expr[gene][ctype] = []
                                for i in idx[ctype]:
                                        expr[gene][ctype].append(int(vals[i]))

        out = open(outname, 'w')
        out.write('Gene\tEndothelial\tOligodendrocytes\tAstrocytes\tMicroglia\tNeurons\tOPC\n')
        for gene in expr:
                outline = [gene]
                for ctype in expr[gene]:
                        c_avg = np.mean(expr[gene][ctype])
                        outline.append(str(c_avg))
                out.write('\t'.join(outline) + '\n')
        out.close()


#Mini mock-test of model build 
#Check to see if cell-type estimates are correct across known fake data
#3 cell types, 3 fake samples, a few fake 'genes' with corresponding fake gene expression that would give known cell type proportions
def bayesTest(mocktable, outname):
	import pymc3 as pymc
        from pymc3.backends import SQLite
        from collections import Counter

        idx = {}
        expr_vector = {}
        for line in open(mocktable):
                if line.startswith('Gene'):
                        header = line.strip().split('\t')
                        for i in range(len(header)):
                                if header[i] != 'Gene':
                                        idx[header[i]] = i
                else:
                        vals = line.strip().split('\t')
                        gene = vals[0]
                        for sample in idx:
                                if sample not in expr_vector:
                                        expr_vector[sample] = [float(vals[idx[sample]])]
                                else:
                                        expr_vector[sample].append(float(vals[idx[sample]]))
        for sample in expr_vector:
                if sample == 'Neurons':
                        neuro = expr_vector[sample]
                if sample == 'Astrocytes':
                        astro = expr_vector[sample]
                if sample == 'Oligodendrocytes':
                        oligo = expr_vector[sample]
                if sample == 'Sample1':
                        one = expr_vector[sample]
                if sample == 'Sample2':
                        two = expr_vector[sample]
                if sample == 'Sample3':
                        three = expr_vector[sample]
        samples = [one, two, three]
        for s in samples:
                model = pymc.Model()
                with pymc.Model() as model:
                        beta = pymc.Dirichlet('beta', a=np.array([1.0, 1.0, 1.0]))
                        sigma = pymc.HalfNormal('sigma', sd=1)
                        y_est = beta[0]*neuro + beta[1]*astro + beta[2]*oligo
                        likelihood = pymc.Normal('y', mu=y_est, sd=sigma, observed=s)
                        trace = pymc.sample(1000, random_seed=123, progressbar=True)
                        s = pymc.summary(trace)
                        #print trace['beta'] #matrix with 3 columns and 1000 rows, need to convert this and do math
                        neurons = trace['beta'][:,0]
                        astrocytes = trace['beta'][:,1]
                        oligodendrocytes = trace['beta'][:,2]
                        n_avg = np.mean(neurons)
                        n_med = np.median(neurons)
                        data = Counter(neurons)
                        data.most_common()
                        n_mode = data.most_common(1)[0][0]
                        print n_avg, n_med, n_mode


#Read in an expression table of genes of interest
#With cell type expression averages as well as sample expression 
#Estimate tissue composition of each sample given the cell type averages across genes of interest
#Genes are rows and samples/cell types are columns
#Header starts with 'Gene'
def CellTypeComp(exprtable, outname):
        import pymc3 as pymc
        from pymc3.backends import SQLite
        from collections import Counter
        import math

        out = open(outname, 'w')

        idx = {}
        samples = {}
        for line in open(exprtable):
                if line.startswith('Gene'):
                        header = line.strip().split('\t')
                        for i in range(len(header)):
                                if header[i] != 'Gene':
                                        idx[header[i]] = i
                else:
                        vals = line.strip().split('\t')
                        if len(vals) == len(header):  #Make sure there is a value for every sample, or don't use
                                for name in idx:
                                        if name not in samples:
                                                samples[name] = [math.log(float(vals[idx[name]])+.1, 10)]
                                        else:
                                                samples[name].append(math.log(float(vals[idx[name]])+.1, 10))

        everything = []
        tot = 0
        array_idx = {}
        #For hg data from Barres lab
        cell_types = ['Endothelial', 'Oligodendrocytes', 'Astrocytes', 'Microglia', 'Neurons', 'OPC']
        for sample in samples:
                if sample not in cell_types:
                        array_idx[sample] = tot
                        everything.append(samples[sample])
                        tot += 1
        neuro = samples['Neurons']
        endo = samples['Endothelial']
        oligo = samples['Oligodendrocytes']
        astro = samples['Astrocytes']
        micro = samples['Microglia']
        opc = samples['OPC']

        everything = np.array(everything)

        for sample in array_idx:
                model = pymc.Model()
                with pymc.Model() as model:
                        beta = pymc.Dirichlet('beta', a=np.array([1.0,1.0,1.0,1.0,1.0]))
                        sigma = pymc.HalfNormal('sigma', sd=1)
                        y_est = beta[0]*neuro + beta[1]*endo + beta[2]*oligo + beta[3]*astro + beta[4]*micro
                        likelihood = pymc.Normal('y', mu=y_est, sd=sigma, observed=everything[array_idx[sample],:])
                        db = SQLite('%s_trace' %(sample))
                        trace = pymc.sample(1000, random_seed=123, progressbar=True, trace=db)
                        s = pymc.summary(trace)
        out.close()


#Wrapper to open all trace files from pymc3 model output in a given directory
#Produce desired stats from each file and output them into a new table
def traceWrapper(exprtable, traceDir, outname):
        import pymc3 as pymc
        from pymc3.backends import SQLite
        from collections import Counter
        import math

        idx = {}
        samples = {}
        for line in open(exprtable):
                if line.startswith('Gene'):
                        header = line.strip().split('\t')
                        for i in range(len(header)):
                                if header[i] != 'Gene':
                                        idx[header[i]] = i
                else:
                        vals = line.strip().split('\t')
                        if len(vals) == len(header):  #Make sure there is a value for every sample, or don't use
                                for name in idx:
                                        if name not in samples:
                                                samples[name] = [math.log(float(vals[idx[name]])+.1, 10)]
                                        else:
                                                samples[name].append(math.log(float(vals[idx[name]])+.1, 10))

        everything = []
        tot = 0
        array_idx = {}
        #For hg:
        cell_types = ['Endothelial', 'Oligodendrocytes', 'Astrocytes', 'Microglia', 'Neurons', 'OPC']
        for sample in samples:
                if sample not in cell_types:
                        array_idx[sample] = tot
                        everything.append(samples[sample])
                        tot += 1
        #For hg:
        neuro = samples['Neurons']
        endo = samples['Endothelial']
        oligo = samples['Oligodendrocytes']
        astro = samples['Astrocytes']
        micro = samples['Microglia']
        opc = samples['OPC']

        everything = np.array(everything)

        out = open(outname, 'w')
        Header = ['Sample', 'Neuronal', 'Endothelial', 'Oligodendrocytes', 'Astrocytes', 'Microglia', 'SE', 'SE_std']
        out.write('\t'.join(Header) + '\n')

        traceDir = os.path.abspath(os.path.expanduser(traceDir))
        files = [f for f in os.listdir(traceDir) if f.endswith('trace')]
        for f in files:

		sample_name = f.strip('_trace')
                path = os.path.join(traceDir, f)

                model = pymc.Model()
                with pymc.Model() as model:
                        beta = pymc.Dirichlet('beta', a=np.array([1.0, 1.0, 1.0, 1.0, 1.0]))
                        sigma = pymc.HalfNormal('sigma', sd=1)
                        y_est = beta[0]*neuro + beta[1]*endo + beta[2]*oligo + beta[3]*astro + beta[4]*micro 
                        likelihood = pymc.Normal('y', mu=y_est, sd=sigma, observed=everything[array_idx[sample],:])
                        trace = pymc.backends.sqlite.load(path)

                outline = [sample_name]

                neurons = trace['beta'][:,0]
                endos = trace['beta'][:,1]
                oligos = trace['beta'][:,2]
                astros = trace['beta'][:,3]
                micros = trace['beta'][:,4]
                types = [neurons, endos, oligos, astros, micros]
                for t in types:
                        avg = np.mean(t)
                        outline.append(avg)
                se = trace['sigma'].mean(axis=0)
                se_std = trace['sigma'].std(axis=0)
                outline.append(se)
                outline.append(se_std)
                out.write('\t'.join(map(str, outline)) + '\n')

        out.close()

