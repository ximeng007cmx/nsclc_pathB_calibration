# References — 全部双源核验通过（2026-10-21）

**核验方式**：`epmc_search.py verify`——① Europe PMC（`EXT_ID:<PMID> AND SRC:MED`，resultType=lite）确认 PMID 真实存在且记录内 DOI 与清单一致；② CrossRef `api.crossref.org/works/<DOI>` 确认 DOI 真实存在。
- 49 条签名文献：`logs/refs_verify_49.log` → **DOI 双源一致 49 条 / 有问题 0 条**（输入 `work/refs_49_to_verify.csv`）。
- 9 条方法学与数据来源文献：**DOI 双源一致 9 条 / 有问题 0 条**（输入 `work/refs_methods_to_verify.csv`）。
- 卷期页码取自 Europe PMC `resultType=core` 记录（`work/refsearch/methods_core.csv`），非凭记忆；期刊缩写按 NLM 惯例。

> ⚠ **错引更正**：Liu 2018（TCGA-CDR）原项目内记录的 DOI `10.1016/j.cell.2018.03.046` 是错的——该 DOI 在 Europe PMC 与 CrossRef 均解析到同刊另一文（*Coupling Neurogenesis to Circuit Formation*, PMID 29625047）。正确 DOI 为 `10.1016/j.cell.2018.02.052`（PMID 29625055、卷期页不变）。pathA 等同源记录如引用该文需同步更正。

## 正文引用（9 条）

1. Yoshihara K, Shahmoradgoli M, Martínez E, Vegesna R, Kim H, Torres-Garcia W, et al. Inferring tumour purity and stromal and immune cell admixture from expression data. Nat Commun;4:2612. PMID 24113773. doi:10.1038/ncomms3612. **[已核验]**
2. Liu J, Lichtenberg T, Hoadley KA, Poisson LM, Lazar AJ, Cherniack AD, et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell;173(2):400-416.e11. PMID 29625055. doi:10.1016/j.cell.2018.02.052. **[已核验]**
3. Barbie DA, Tamayo P, Boehm JS, Kim SY, Moody SE, Dunn IF, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature;462(7269):108-112. PMID 19847166. doi:10.1038/nature08460. **[已核验]**
4. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics;14:7. PMID 23323831. doi:10.1186/1471-2105-14-7. **[已核验]**
5. Becht E, Giraldo NA, Lacroix L, Buttard B, Elarouci N, Petitprez F, et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biol;17(1):218. PMID 27765066. doi:10.1186/s13059-016-1070-5. **[已核验]**
6. Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biol;18(1):220. PMID 29141660. doi:10.1186/s13059-017-1349-1. **[已核验]**
7. Director's Challenge Consortium for the Molecular Classification of Lung Adenocarcinoma, Shedden K, Taylor JM, Enkemann SA, Tsao MS, Yeatman TJ, et al. Gene expression-based survival prediction in lung adenocarcinoma: a multi-site, blinded validation study. Nat Med;14(8):822-827. PMID 18641660. doi:10.1038/nm.1790. **[已核验]**
8. Okayama H, Kohno T, Ishii Y, Shimada Y, Shiraishi K, Iwakawa R, et al. Identification of genes upregulated in ALK-positive and EGFR/KRAS/ALK-negative lung adenocarcinomas. Cancer Res;72(1):100-111. PMID 22080568. doi:10.1158/0008-5472.can-11-1403. **[已核验]**
9. Athar A, Füllgrabe A, George N, Iqbal H, Huerta L, Ali A, et al. ArrayExpress update - from bulk to single-cell expression data. Nucleic Acids Res;47(D1):D711-D715. PMID 30357387. doi:10.1093/nar/gky964. **[已核验]**

## 签名文献 49 条（Table 1，逐条 PMID+DOI 核验）

- **SIG01** [inflammation] An Inflammation-Related Nine-Gene Signature to Improve Prognosis Prediction of Lung Adenocarcinoma. Dis Markers;2021. PMID 34580602. doi:10.1155/2021/9568057. **[已核验]**
- **SIG02** [inflammation] A Novel Inflammatory-Related Gene Signature Based Model for Risk Stratification and Prognosis Prediction in Lung Adenocarcinoma. Front Genet;2021. PMID 35069695. doi:10.3389/fgene.2021.798131. **[已核验]**
- **SIG03** [inflammation] Development and Validation of a 7-Gene Inflammatory Signature Forecasts Prognosis and Diverse Immune Landscape in Lung Adenocarcinoma. Front Mol Biosci;2022. PMID 35372503. doi:10.3389/fmolb.2022.822739. **[已核验]**
- **SIG04** [inflammation] Identification of an inflammatory response signature associated with prognostic stratification and drug sensitivity in lung adenocarcinoma. Sci Rep;2022. PMID 35710585. doi:10.1038/s41598-022-14323-6. **[已核验]**
- **SIG05** [inflammation] Identification and validation of a novel prognostic model of inflammation-related gene signature of lung adenocarcinoma. Sci Rep;2022. PMID 36042374. doi:10.1038/s41598-022-19105-8. **[已核验]**
- **SIG06** [inflammation] Risk stratification and prognosis prediction based on inflammation-related gene signature in lung squamous carcinoma. Cancer Med;2023. PMID 36056909. doi:10.1002/cam4.5190. **[已核验]**
- **SIG07** [inflammation] Two inflammation-related genes model could predict risk in prognosis of patients with lung adenocarcinoma. Clin Transl Oncol;2025. PMID 40053252. doi:10.1007/s12094-025-03861-w. **[已核验]**
- **SIG08** [inflammation] Identification of Inflammatory Cell Death-Related Prognostic Signature and Its Regulatory Mechanism in Lung Adenocarcinoma. Cancer Med;2026. PMID 42760759. doi:10.1002/cam4.72255. **[已核验]**
- **SIG09** [metabolism] A Distinct Glucose Metabolism Signature of Lung Adenocarcinoma With Prognostic Value. Front Genet;2022. PMID 35615380. doi:10.3389/fgene.2022.860677. **[已核验]**
- **SIG10** [metabolism] Fatty acid metabolism prognostic signature predicts tumor immune microenvironment and immunotherapy, and identifies tumorigenic role of MOGAT2 in lung adenocarcinoma. Front Immunol;2024. PMID 39478862. doi:10.3389/fimmu.2024.1456719. **[已核验]**
- **SIG11** [metabolism] Clinical Significance and Immunometabolism Landscapes of a Novel Recurrence-Associated Lipid Metabolism Signature In Early-Stage Lung Adenocarcinoma: A Comprehensive Analysis. Front Immunol;2022. PMID 35222371. doi:10.3389/fimmu.2022.783495. **[已核验]**
- **SIG12** [metabolism] Transcriptome-Based Six-Gene Fatty Acid Metabolism Signature for Prognosis and Predicted Immunotherapy Response in Lung Adenocarcinoma: Cross-Population Validation. Genes (Basel);2026. PMID 42650107. doi:10.3390/genes17080914. **[已核验]**
- **SIG13** [metabolism] Integrating bulk and single-cell transcriptomic data to construct a risk model for histidine metabolism-related epithelial cell features in lung adenocarcinoma, predicting prognosis and immune landscape. RNA Biol;2026. PMID 42003422. doi:10.1080/15476286.2026.2662721. **[已核验]**
- **SIG14** [metabolism] Prognostic Implication of Energy Metabolism-Related Gene Signatures in Lung Adenocarcinoma. Front Oncol;2022. PMID 35494074. doi:10.3389/fonc.2022.867470. **[已核验]**
- **SIG15** [metabolism] Development of a gene signature associated with iron metabolism in lung adenocarcinoma. Bioengineered;2021. PMID 34323652. doi:10.1080/21655979.2021.1954840. **[已核验]**
- **SIG16** [metabolism] Lipid metabolism-related gene expression predicts prognostic outcomes in lung adenocarcinoma. Transl Cancer Res;2026. PMID 42180928. doi:10.21037/tcr-2025-1-2665. **[已核验]**
- **SIG17** [metabolism] Construction and validation of an oxidative phosphorylation-related gene signature in lung squamous cell carcinoma patients. Lung Cancer Manag;2025. PMID 41343408. doi:10.1080/17581966.2025.2597179. **[已核验]**
- **SIG18** [metabolism] Establishment and validation of an eight-gene metabolic-related prognostic signature model for lung adenocarcinoma. Aging (Albany NY);2021. PMID 33619235. doi:10.18632/aging.202681. **[已核验]**
- **SIG19** [metabolism] Development and internal validation of a six-gene prognostic model based on galactose metabolism for overall survival in lung adenocarcinoma. Transl Cancer Res;2026. PMID 42724461. doi:10.21037/tcr-2026-0696. **[已核验]**
- **SIG20** [metabolism] Lipid metabolism gene-wide profile and survival signature of lung adenocarcinoma. Lipids Health Dis;2020. PMID 33050938. doi:10.1186/s12944-020-01390-9. **[已核验]**
- **SIG21** [metabolism] Identification of a six-gene metabolic signature predicting overall survival for patients with lung adenocarcinoma. PeerJ;2020. PMID 33344071. doi:10.7717/peerj.10320. **[已核验]**
- **SIG22** [metabolism] Identification of a Prognostic Gene Signature for Chemoresistance Prediction in Lung Adenocarcinoma by Screening Mitochondrial Metabolism Gene Sets. Int J Mol Sci;2026. PMID 41977253. doi:10.3390/ijms27073065. **[已核验]**
- **SIG23** [metabolism] Multiomic and functional validation of ACSL3, a regulator of fatty acid metabolism, as a lymph node metastasis-associated gene in lung adenocarcinoma. Respir Res;2026. PMID 41998631. doi:10.1186/s12931-026-03671-w. **[已核验]**
- **SIG24** [metabolism] Integrative analysis of a novel signature incorporating metabolism and stemness-related genes for risk stratification and assessing clinical outcomes and therapeutic responses in lung adenocarcinoma. BMC Cancer;2025. PMID 40170009. doi:10.1186/s12885-025-13984-6. **[已核验]**
- **SIG25** [metabolism] Predicting Differences in Treatment Response and Survival Time of Lung Adenocarcinoma Patients Based on a Prognostic Risk Model of Glycolysis-Related Genes. Front Genet;2022. PMID 35692818. doi:10.3389/fgene.2022.828543. **[已核验]**
- **SIG26** [metabolism] Multi-omics analysis and experiments validate the tumor-suppressive role of mitochondrial lipid metabolism gene ACSM5 in lung adenocarcinoma and its impact on the immune microenvironment. Hum Genomics;2026. PMID 41634895. doi:10.1186/s40246-026-00916-z. **[已核验]**
- **SIG27** [metabolism] ‌PPIA regulates fatty acid and glutamine metabolism in lung adenocarcinoma based on multiomics prognostic model and experiment validation. Sci Rep;2026. PMID 41491011. doi:10.1038/s41598-025-34313-8. **[已核验]**
- **SIG28** [metabolism] Applying machine learning algorithms to develop a survival prediction model for lung adenocarcinoma based on genes related to fatty acid metabolism. Front Pharmacol;2023. PMID 37920207. doi:10.3389/fphar.2023.1260742. **[已核验]**
- **SIG29** [metabolism] Characteristics of folic acid metabolism-related genes unveil prognosis and treatment strategy in lung adenocarcinoma. BMC Pulm Med;2025. PMID 40405133. doi:10.1186/s12890-025-03694-x. **[已核验]**
- **SIG30** [metabolism] Construction of lung adenocarcinoma subtype and prognosis model based on fatty acid metabolism-related genes. Discov Oncol;2025. PMID 40405049. doi:10.1007/s12672-025-02613-0. **[已核验]**
- **SIG31** [metabolism] An abnormal metabolism-related gene, ALG3, is a potential diagnostic and prognostic biomarker for lung adenocarcinoma. Medicine (Baltimore);2024. PMID 39287231. doi:10.1097/md.0000000000038746. **[已核验]**
- **SIG32** [inflammation] Clinical Significance and Immunologic Landscape of a Five-IL(R)-Based Signature in Lung Adenocarcinoma. Front Immunol;2021. PMID 34497605. doi:10.3389/fimmu.2021.693062. **[已核验]**
- **SIG33** [metabolism] Choline metabolism reprogramming mediates an immunosuppressive microenvironment in non-small cell lung cancer (NSCLC) by promoting tumor-associated macrophage functional polarization and endothelial cell proliferation. J Transl Med;2024. PMID 38730286. doi:10.1186/s12967-024-05242-3. **[已核验]**
- **SIG34** [metabolism] A Novel Integrated Metabolism-Immunity Gene Expression Model Predicts the Prognosis of Lung Adenocarcinoma Patients. Front Pharmacol;2021. PMID 34393804. doi:10.3389/fphar.2021.728368. **[已核验]**
- **SIG35** [metabolism] The systematic analysis of genes related to butyrate metabolism suggests that CDKN3 could serve as a promising therapeutic target for lung adenocarcinoma treatment. BMC Cancer;2025. PMID 39806313. doi:10.1186/s12885-024-13409-w. **[已核验]**
- **SIG36** [metabolism] A new prognostic model for lung adenocarcinoma according propionate metabolism related genes: a comprehensive bioinformatic study. Discov Oncol;2025. PMID 41003841. doi:10.1007/s12672-025-03540-w. **[已核验]**
- **SIG37** [metabolism] Insight of a Metabolic Prognostic Model to Identify Tumor Environment and Drug Vulnerability for Lung Adenocarcinoma. Front Immunol;2022. PMID 35812404. doi:10.3389/fimmu.2022.872910. **[已核验]**
- **SIG38** [metabolism] Prognostic implication of glycolysis related gene signature in non-small cell lung cancer. J Cancer;2021. PMID 33403045. doi:10.7150/jca.50274. **[已核验]**
- **SIG39** [metabolism] Identification of a 5-Gene Metabolic Signature for Predicting Prognosis Based on an Integrated Analysis of Tumor Microenvironment in Lung Adenocarcinoma. J Oncol;2020. PMID 32684932. doi:10.1155/2020/5310793. **[已核验]**
- **SIG40** [metabolism] ANLN: A New Hub in Glutamine Metabolism of Lung Adenocarcinoma by scRNA-Seq and Machine Learning. Oncol Res;2026. PMID 42500572. doi:10.32604/or.2026.079515. **[已核验]**
- **SIG41** [metabolism] Establishment of the Prognostic Index Reflecting Tumor Immune Microenvironment of Lung Adenocarcinoma Based on Metabolism-Related Genes. J Cancer;2020. PMID 33193873. doi:10.7150/jca.49266. **[已核验]**
- **SIG42** [inflammation] Identifying Chemokine System-Related Phenotype to Predict Immune Feature in Pan-Cancer and Prognostic Signature for Lung Adenocarcinoma. J Inflamm Res;2025. PMID 41017970. doi:10.2147/jir.s537256. **[已核验]**
- **SIG43** [inflammation] A novel gene signature derived from the CXC subfamily of chemokine receptors predicts the prognosis and immune infiltration of patients with lung adenocarcinoma. Medicine (Baltimore);2022. PMID 36254009. doi:10.1097/md.0000000000030982. **[已核验]**
- **SIG44** [immune] Identification of immune-related prognostic biomarkers in lung squamous cell carcinoma microenvironment. Front Immunol;2025. PMID 41562084. doi:10.3389/fimmu.2025.1724319. **[已核验]**
- **SIG45** [immune] Identification and Validation of Immune-Related Gene Signature for Predicting Lymph Node Metastasis and Prognosis in Lung Adenocarcinoma. Front Mol Biosci;2021. PMID 34109216. doi:10.3389/fmolb.2021.679031. **[已核验]**
- **SIG46** [immune] Identification of an Innate Immune-Related Prognostic Signature in Early-Stage Lung Squamous Cell Carcinoma. Int J Gen Med;2021. PMID 34876838. doi:10.2147/ijgm.s341175. **[已核验]**
- **SIG47** [immune] Development of an Immune-Related Risk Signature for Predicting Prognosis in Lung Squamous Cell Carcinoma. Front Genet;2020. PMID 33005178. doi:10.3389/fgene.2020.00978. **[已核验]**
- **SIG48** [immune] Immune-Related RNA-Binding Protein-Based Signature With Predictive and Prognostic Implications in Patients With Lung Adenocarcinoma. Front Mol Biosci;2022. PMID 35647031. doi:10.3389/fmolb.2022.807622. **[已核验]**
- **SIG49** [immune] Identification of a costimulatory molecule-based signature for predicting prognosis risk and immunotherapy response in patients with lung adenocarcinoma. Oncoimmunology;2020. PMID 33457102. doi:10.1080/2162402x.2020.1824641. **[已核验]**
