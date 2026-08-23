"""Fixture copies of the VOGDB tables, with the real headers and row shapes.

Kept in one place so the parser tests and the source-integration test share them
and neither needs the 554 MB profile archive (see ``proposal.md`` §10.1).
"""

ANNOTATIONS = """\
#GroupName\tProteinCount\tSpeciesCount\tFunctionalCategory\tConsensusFunctionalDescription
VOG00001\t10\t10\tXu\tREFSEQ hypothetical protein
VOG00002\t8\t7\tXr\tglycoprotein
VOG00003\t2\t2\tXs\tcapsid protein
VOG00004\t287\t285\tXh\tterminase large subunit
VOG00005\t12\t11\tXr\tRNA-dependent RNA polymerase
"""

LCA = """\
#GroupName\tGenomesInGroupAndLCA\tGenomesTotalInLCA\tLastCommonAncestor_Name\tLastCommonAncestor_Taxon_ID
VOG00001\t10\t14\tViruses;Riboviria;Chuviridae;Mivirus\t111111
VOG00002\t8\t8\tViruses;Riboviria;Chuviridae;Mivirus\t222222
VOG00003\t2\t2\tViruses;Riboviria;Flaviviridae\t333333
VOG00004\t285\t5444\tViruses;Duplodnaviria;Caudoviricetes\t2731619
VOG00005\t12\t12\tViruses;Riboviria;Chuviridae;Culicidavirus\t444444
"""

MEMBERS = """\
#GroupName\tProteinCount\tSpeciesCount\tFunctionalCategory\tProteinIDs
VOG00001\t10\t10\tXu\t1001.YP_000001.1,1002.YP_000002.1
VOG00002\t8\t7\tXr\t1001.YP_000003.1,1002.YP_000004.1
VOG00003\t2\t2\tXs\t1003.YP_000005.1
VOG00004\t287\t285\tXh\t2001.YP_000006.1,2002.YP_000007.1
VOG00005\t12\t11\tXr\t1001.YP_000008.1,1004.YP_000009.1
"""

VIRUSONLY = """\
#GroupName\tOnly in viruses (high stringency)\tOnly in viruses (medium stringency)\tOnly in viruses (low stringency)
VOG00001\t1\t1\t1
VOG00002\t1\t1\t1
VOG00003\t1\t1\t1
VOG00004\t0\t0\t0
VOG00005\t0\t1\t1
"""

SPECIES = """\
#species name\ttaxon id\tsource\tsource version
Mivirus chuvi\t1001\tNCBI Refseq\t236
Mivirus alpha\t1002\tNCBI Refseq\t236
Zika virus\t1003\tNCBI Refseq\t236
Escherichia phage T4\t2001\tNCBI Refseq\t236
Escherichia phage T7\t2002\tNCBI Refseq\t236
Culicidavirus imjinense\t1004\tNCBI Refseq\t236
"""

HOSTS = """\
#taxon id\tphage/nonphage\thost\tsuperkingdom of host
1001\tnonphage\tAedes aegypti\t
1002\tnonphage\tAedes aegypti\t
1003\tnonphage\tHomo sapiens\t
1004\tnonphage\t\t
2001\tphage\t\t
2002\tphage\t\t
"""


#: ``filename -> table`` in the layout a VOGDB download uses.
VOGDB_FIXTURE_TABLES = {
    "vog.annotations.tsv.gz": ANNOTATIONS,
    "vog.lca.tsv.gz": LCA,
    "vog.members.tsv.gz": MEMBERS,
    "vog.virusonly.tsv.gz": VIRUSONLY,
    "vogdb.species.txt": SPECIES,
    "vogdb.host.txt": HOSTS,
}
