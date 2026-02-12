import pandas as pd
import cobra
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"
import networkx as nx
from collections import defaultdict
from collections import Counter
import pathlib



out_dir = "tutorial_output/construct_recon3d_met_net/"
pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

# Read genome scale metabolic model
met_model = cobra.io.load_matlab_model("tutorial_data/Recon3D_301.mat")


met_df = pd.DataFrame(
    [(x.id, x.name) for x in met_model.metabolites],
    columns=["metabolite_id", "metabolite_name"])

met_df.to_csv(out_dir + "recon3d_metabolites.tsv", sep="\t", index=False)


# If the genome scale metabolic model does not have reaction.subsystem
# attribute, check if the pathway information is stored in other reaction
# attributes.
rxn_gene_df = []
for i_rxn in met_model.reactions:
    for j in i_rxn.genes:
        rxn_gene_df.append(
            (i_rxn.id, i_rxn.name, j.id, str(i_rxn), i_rxn.subsystem))


# Delete iterator variables to avoid using them outside of the loops by mistake.
del i_rxn, j

rxn_gene_df = pd.DataFrame(
    sorted(rxn_gene_df),
    columns=["rxn", "rxn_name", "gene_id", "equation", "pathway"])



# Count the number of metabolites for each compartment.
mets = []
compartment_counter = []
for i in met_model.reactions:
    for j in i.metabolites:
        mets.append(j.id)
        compartment_counter.append("[" + j.compartment + "]")


del i, j

mets = Counter(mets)
compartment_counter = Counter(compartment_counter)



# Manually selected pervasive metabolites.
f_mets = set([
    "h2o", "atp", "pi", "adp", "na1", "nadp", "nadph", "coa", "o2", "hco3",
    "crn", "ppi", "nad", "amp", "nadh", "co2", "udp", "h2o2", "nh4", "cl",
    "accoa", "ppi", "ppcoa", "fad", "fadh2", "occoa", "amp", "chsterol", "so4",
])


# Create a dictionary of edges and weights.
edge_lut = defaultdict(int)

for i_rxn in met_model.reactions:
    i_mets = set([x.id for x in i_rxn.metabolites.keys()])
    for j_rxn in met_model.reactions:
        j_mets = [x.id for x in j_rxn.metabolites.keys()]
        weight = 0
        for k in j_mets:
            if (k in i_mets) and (k[:-3] not in f_mets):
                weight += 1 / (mets[k] / compartment_counter[k[-3:]])
        edge_key = tuple(sorted([i_rxn.id, j_rxn.id]))
        if edge_key not in edge_lut:
            edge_lut[edge_key] = weight
        else:
            assert weight == edge_lut[edge_key]


del i_rxn, i_mets, j_rxn, j_mets, weight, k, edge_key


edge_df = pd.DataFrame(
    [(k[0], k[1], v) for k, v in edge_lut.items()],
    columns=["src", "dest", "mn_weight"]
)


nzw_edge_df = edge_df.loc[edge_df["mn_weight"] != 0, ]


# Map Entrez Gene ID to gene symbols.
#
# The mappings were obtained using BioMart.
eid_gs_df = pd.read_table(
    "tutorial_data/biomart_ensembl_110_eid_hs_mapping.tsv",
    dtype=str)

eid_gs_dict = dict(zip(eid_gs_df.gene_entrez_id, eid_gs_df.gene))


gene_symbols = []
for i in rxn_gene_df.gene_id.tolist():
    i_dot_pos = i.find(".")
    if  i_dot_pos >= 0:
        i_eid = i[:i_dot_pos]
    else:
        i_eid = i
    gene_symbols.append(eid_gs_dict.get(i_eid, ""))


del i, i_dot_pos, i_eid

m_rxn_gene_df = rxn_gene_df.copy()
m_rxn_gene_df["gene"] = gene_symbols

m_rxn_gene_df = m_rxn_gene_df.loc[
    ~m_rxn_gene_df.gene.str.fullmatch(""), :]

m_rxn_gene_df = m_rxn_gene_df.loc[
    :, ["rxn", "rxn_name", "equation", "pathway", "gene"]]

common_rxns = set.intersection(
    set(m_rxn_gene_df.rxn.tolist()),
    set(nzw_edge_df.src.tolist()),
    set(nzw_edge_df.dest.tolist())
)


def drop_dup_sort_rows(df):
    df = df.copy()
    df = df.drop_duplicates()
    df = df.sort_values(by=df.columns.tolist())
    return df


output_rxn_gene_df = drop_dup_sort_rows(
    m_rxn_gene_df.loc[m_rxn_gene_df.rxn.isin(common_rxns), :])

output_nzw_edge_df = drop_dup_sort_rows(
    nzw_edge_df.loc[
        nzw_edge_df.src.isin(common_rxns) &
            nzw_edge_df.dest.isin(common_rxns), :]
)

output_rxn_gene_df.to_csv(
    out_dir + "recon3d_rxn_gene.tsv", sep="\t", index=False)

output_nzw_edge_df.to_csv(
    out_dir + "recon3d_edges.tsv.gz", sep="\t", index=False,
    compression={"method": "gzip", "compresslevel": 1, "mtime": 1})
