# Literature search behind Section 1.3 of the Physiological Measurement manuscript

This file deposits the search strings and the screening lists referred to in Section 1.3 and in the Data availability statement of
*Measurement-domain component removal in thoracic electrical impedance tomography reconstructs the removed component, not the heart: a phantom study with a human recording* (Lin and Ko, submitted to *Physiological Measurement*, 2026).

The searches were run on **24 August 2026**. Their queries and counts were recorded in the manuscript text at the time; the per-item screening lists were not saved as files then. They were therefore **re-run on 31 October 2026** with the same queries, and every item returned today is listed below with its screening decision. Items that appeared after 24 August 2026 are marked. The machine-readable listing is `literature_search_rerun_2026-10-31.json`.

Both statements that rest on this search are negative and dated: that DAGAF had not been applied to EIT, and that no two measurement-domain decompositions had been compared under one contract, in the sources searched, as of 24 August 2026. Neither Scopus nor Web of Science was searched.

## 1. PubMed

Query 1 — DAGAF applied to EIT:

```
("electrical impedance tomography"[tiab] OR EIT[tiab]) AND ("Gaussian average filtering decomposition"[tiab] OR DAGAF[tiab] OR GAFD[tiab])
```

24 August 2026: **0 records**. 31 October 2026: **0 records**.

Query 2 — decomposition for cardiac or perfusion signals in EIT:

```
"electrical impedance tomography"[tiab] AND decomposition[tiab] AND (cardiac[tiab] OR perfusion[tiab] OR heart[tiab])
```

24 August 2026: **9 records**. 31 October 2026: **10 records**; the addition is PMID 42710529 (Hidayah et al 2026, published 17 September 2026). Screening of the ten:

| PMID | Date | Journal | Title | Decision |
|---|---|---|---|---|
| 42710529 | 2026 Sep 17 | Biomedical physics & engineering express | Cardiac and respiratory signal separation in electrical impedance tomography (EIT). | not DAGAF; after the search date; see the OpenAlex note on the same paper below |
| 40288963 | 2025 Apr 25 | Sheng wu yi xue gong cheng xue za zhi = Journal of biomedical engineering = Shengwu yixue gongchengxue zazhi | [Study on the separation method of lung ventilation and lung perfusion signals in electrical impedance tomography based on rime algorithm optimized variational mode decomposition]. | not DAGAF; VMD, one method (journal version of the work also reported at a conference) |
| 38925138 | 2024 Jul 26 | Physiological measurement | Separation of ventilation and perfusion of electrical impedance tomography image streams using multi-dimensional ensemble empirical mode decomposition. | not DAGAF; multi-dimensional ensemble EMD on image streams; one method, image domain |
| 38697210 | 2024 May 21 | Physiological measurement | Improved filtering methods to suppress cardiovascular contamination in electrical impedance tomography recordings. | not DAGAF; filtering of cardiovascular contamination; no decomposition comparison |
| 35570528 | 2022 | Current medical imaging | Separation of Heart and Lung-related Signals in Electrical Impedance Tomography Using Empirical Mode Decomposition. | not DAGAF; the previous study itself (reference 19 of the manuscript) |
| 35174192 | 2022 | Frontiers in medicine | Removing Clinical Motion Artifacts During Ventilation Monitoring With Electrical Impedance Tomography: Introduction of Methodology and Validation With Simulation and Patient Data. | not DAGAF; motion-artefact removal; not cardiac separation |
| 34167091 | 2021 Aug 27 | Physiological measurement | A model-based source separation algorithm for lung perfusion imaging using electrical impedance tomography. | not DAGAF; model-based source separation; one method |
| 33529234 | 2021 | PloS one | Semi-Siamese U-Net for separation of lung and heart bioimpedance images: A simulation study of thorax EIT. | not DAGAF; image-domain deep learning (reference 18 of the manuscript) |
| 33021240 | 2020 Oct 6 | Physiological measurement | Cardiac eigen imaging: a novel method to isolate cardiac activity in thoracic electrical impedance tomography. | not DAGAF; cardiac eigen imaging, PCA on the image series; one method, image domain |
| 9001603 | 1996 Nov | Physiological measurement | An interleaved drive electrical impedance tomography image reconstruction algorithm. | not DAGAF; 1996 reconstruction algorithm; not relevant |

None of the ten applies DAGAF to EIT, and none compares two single-frame measurement-domain decompositions under one additive contract.

## 2. OpenAlex forward citations

### 2.1 Works citing the DAGAF paper (Lin, Tan and Tian 2022, *Biomed Signal Process Control* 71:103104; reference 21)

OpenAlex id `W3198757505`. 24 August 2026: **29 citing works**. 31 October 2026: **31**.

Screening question: is the citing work an application of DAGAF to electrical impedance tomography? Decision for all 31: **no**. None is in EIT; the nearest is an impedance cardiography paper (Ma et al 2026, *Measurement*), a single-channel modality, which was already among the 29 in August. One of the two additions since August carries a publication date after 24 August (Indrasree et al, 19 September 2026, lupus nephritis image analysis); the other cannot be identified from dates alone, because OpenAlex also indexes works with earlier dates late.

| # | Date | Type | Venue | First author | Title | DOI | Decision |
|---|---|---|---|---|---|---|---|
| 1 | 2022-11-07 | conference-paper | Asia-Pacific Signal and Information Processing Association Annual Summit and Conference (APSIPA ASC) | Jin Wu et al | Design of Optimal FIR Digital Filter by Swarm Optimization Technique | 10.23919/apsipaasc55919.2022.9980121 | exclude: not EIT |
| 2 | 2023-01-01 | article | Journal of Instrumentation | Ziwei Tian et al | A method for portable PPG signal dynamic acquisition and denoise based on CEEMDAN and wavelet transform | 10.1088/1748-0221/18/01/p01023 | exclude: not EIT |
| 3 | 2023-01-01 | article | Chinese Journal of Lasers | 陶嘉庆 Tao Jiaqing et al | 基于多区域主频增强方法获取颈部的弱脉搏信号 | 10.3788/cjl221273 | exclude: not EIT |
| 4 | 2023-01-01 | conference-paper | Advances in intelligent systems and computing | Maksim Vladimirovich Bobyr et al | A Two-Level Fuzzy Model for Filtering Signals of the Automatic Dependent Surveillance-Broadcast | 10.1007/978-981-99-1472-2_5 | exclude: not EIT |
| 5 | 2023-04-07 | article | Sensors | Yue‐Der Lin et al | A Frequency Estimation Scheme Based on Gaussian Average Filtering Decomposition and Hilbert Transform: With Estimation of Respiratory Rate as an Example | 10.3390/s23083785 | exclude: not EIT |
| 6 | 2023-05-26 | article | Journal of Nondestructive Evaluation | Clément Remacha et al | Eigen Background Subtraction for Industrial Flaw Detection: Application to High-Pressure Turbine Blade CT Scans | 10.1007/s10921-023-00955-9 | exclude: not EIT |
| 7 | 2023-10-06 | article | IEEE Sensors Journal | Kuan-Chung Ting et al | Detection of Otitis Media With Effusion Using In-Ear Microphones and Machine Learning | 10.1109/jsen.2023.3321093 | exclude: not EIT |
| 8 | 2023-12-19 | article | Proceedings of Southwest State University | A. E. Arkhipov | Filtering of Complex Signals Based on a Two-Level Fuzzy-Logic Model | 10.21869/2223-1560-2023-27-2-140-154 | exclude: not EIT |
| 9 | 2024-04-04 | article | Journal of Electrical Systems | Xiaohua Su | English Teaching Quality Evaluation Model for Higher Education Based on Educational Data Mining Technology | 10.52783/jes.3150 | exclude: not EIT |
| 10 | 2024-04-04 | article | Journal of Electrical Systems | Xinye Wang | Augmented Physics-Informed Neural Networks and Seahorse Optimization Algorithm for Music Scene Simulation and Emotion Regulation in Virtual Reality Technology | 10.52783/jes.3060 | exclude: not EIT |
| 11 | 2024-04-04 | article | Journal of Electrical Systems | Ruifen Zhao | Construction and Optimization of Financial Transformation Decision Model Based on Data Mining | 10.52783/jes.3142 | exclude: not EIT |
| 12 | 2024-04-04 | article | Journal of Electrical Systems | Jinqing Luo | Optimization Research on Provincial Users' Participation in Regional Electricity Spot Market Bidding Based on Hierarchically Gated Recurrent Neural Network | 10.52783/jes.3167 | exclude: not EIT |
| 13 | 2024-05-11 | article | Remote Sensing | Baofeng Tian et al | Envelope Extraction Algorithm for Magnetic Resonance Sounding Signals Based on Adaptive Gaussian Filters | 10.3390/rs16101713 | exclude: not EIT |
| 14 | 2024-06-14 | article | Photonics | Ningning Han et al | Modulated Short-Time Fourier-Transform-Based Nonstationary Signal Decomposition for Dual-Comb Ranging Systems | 10.3390/photonics11060560 | exclude: not EIT |
| 15 | 2024-08-08 | article | Knowledge-Based Systems | C. Moorthy et al | Optimized brain tumor identification via graph sample and aggregate-attention network with Artificial Lizard Search Algorithm | 10.1016/j.knosys.2024.112362 | exclude: not EIT |
| 16 | 2024-09-06 | article | Knowledge-Based Systems | Atul B. Kathole et al | Solar energy prediction in IoT system based optimized complex-valued spatio-temporal graph convolutional neural network | 10.1016/j.knosys.2024.112400 | exclude: not EIT |
| 17 | 2024-10-03 | article | Information | Yue‐Der Lin et al | Fall Detection Based on Data-Adaptive Gaussian Average Filtering Decomposition and Machine Learning | 10.3390/info15100606 | exclude: not EIT |
| 18 | 2024-10-18 | article | Biomedical Signal Processing and Control | S. Perumal Sankar et al | Optimized neural network for vulnerable plaque detection in OCT images with noise tolerance and adaptive coefficient zeroing | 10.1016/j.bspc.2024.107046 | exclude: not EIT |
| 19 | 2024-10-30 | article | IEEE Sensors Journal | Kai-Chun Liu et al | A Cross-Modal Autoencoder for Contactless Electrocardiography Monitoring Using Frequency-Modulated Continuous Wave Radar | 10.1109/jsen.2024.3486154 | exclude: not EIT |
| 20 | 2024-12-14 | article | Biomedical Signal Processing and Control | Monika A. Prucnal et al | Improving the quality of pulse rate variability derived from wearable devices using adaptive, spectrum and nonlinear filtering | 10.1016/j.bspc.2024.107336 | exclude: not EIT |
| 21 | 2024-12-25 | article | Pattern Recognition | Moorthy V et al | Deep vision-based wildlife intrusion detection with colour distribution preserved generative adversarial networks | 10.1016/j.patcog.2024.111272 | exclude: not EIT |
| 22 | 2025-05-09 | article | Biomedical Signal Processing and Control | Tien-Hsiung Ku et al | A computationally efficient algorithm for estimating respiratory rate from seismocardiogram | 10.1016/j.bspc.2025.108030 | exclude: not EIT |
| 23 | 2025-05-09 | conference-paper | International Conference on Energy, Power and Environment (ICEPE) | G. Kalyana Chakravarthy et al | An Advanced EEG based Human Emotion Detection with Giant Armadillo Optimization Algorithm | 10.1109/icepe65965.2025.11139644 | exclude: not EIT |
| 24 | 2025-08-28 | article | Neuromorphic Computing and Engineering | Yuqi Ding et al | Physical reservoir computing for biomedical applications | 10.1088/2634-4386/ae006c | exclude: not EIT |
| 25 | 2025-10-24 | conference-paper | International Conference on Communication, Computer, and Information Technology (IC3IT) | G. Deepika et al | Adaptive Lightweight Autoencoder with Noise Estimation Module for Noise Reduction in ECG Signals | 10.1109/ic3it66137.2025.11340876 | exclude: not EIT |
| 26 | 2025-11-20 | article | Pattern Recognition Letters | Samuele Salti et al | Additive decomposition of one-dimensional signals using Transformers | 10.1016/j.patrec.2025.11.002 | exclude: not EIT |
| 27 | 2025-12-01 | article | Biomedical Signal Processing and Control | Kai-Chun Liu et al | Efficient fall detection using Kalman filter-enhanced triaxial accelerometer signals and machine learning | 10.1016/j.bspc.2025.109304 | exclude: not EIT |
| 28 | 2026-01-01 | book-chapter | Elsevier eBooks | Ildar Rakhmatulin et al | Signal processing for biomedical applications | 10.1016/b978-0-443-49276-1.00015-3 | exclude: not EIT |
| 29 | 2026-02-01 | article | Automation and Remote Control | D. Viji et al | Droop Control-Based Optimization of Hybrid Renewable Energy Systems: Combining Node-Level Capsule Graph Neural Network with Artificial Protozoa Optimizer | 10.1134/s0005117925601204 | exclude: not EIT |
| 30 | 2026-06-08 | article | Measurement | Shuai MA et al | Lightweight attention diffusion model for denoising and reconstruction of impedance cardiography signals | 10.1016/j.measurement.2026.122202 | exclude: not EIT; nearest (impedance cardiography) |
| 31 | 2026-09-19 **(after 24 Aug 2026)** | article | International Urology and Nephrology | Gogineni Indrasree et al | An explainable hybrid SMPR-Net-XGBoost framework for automated lupus nephritis detection using medical image analysis | 10.1007/s11255-026-05386-9 | exclude: not EIT |

### 2.2 Works citing the previous measurement-domain study (Cheng, Su and Ko 2022, *Curr Med Imaging* 18:1396; reference 19)

OpenAlex id `W4280513786`. 24 August 2026: **9 citing works**. 31 October 2026: **10**.

Screening question: does the citing work compare two or more measurement-domain decompositions of a single EIT frame under one contract? Decision for all ten: **no**. Exactly one work was added since August, Hidayah et al 2026 (8 September 2026), which is the same paper that PubMed query 2 now returns; its note is below.

| # | Date | Type | Venue | First author | Title | DOI | Decision |
|---|---|---|---|---|---|---|---|
| 1 | 2022-10-20 | article | Physiological Measurement | Ke Zhang et al | Deep feature-domain matching for cardiac-related component separation from a chest electrical impedance tomography image series: proof-of-concept study | 10.1088/1361-6579/ac9c44 | exclude: EIT, cardiac-related component; unsupervised deep feature-domain matching on image series; one method, image domain |
| 2 | 2024-05-13 | review | La radiologia medica | Honghao Song et al | The influence of manual segmentation strategies and different phases selection on machine learning-based computed tomography in renal tumors: a systematic review and meta-analysis | 10.1007/s11547-024-01825-8 | exclude: renal CT review; not EIT |
| 3 | 2024-06-26 | article | Physiological Measurement | Alfred Hülkenberg et al | Separation of ventilation and perfusion of electrical impedance tomography image streams using multi-dimensional ensemble empirical mode decomposition | 10.1088/1361-6579/ad5c39 | exclude: EIT; multi-dimensional ensemble EMD on image streams; one method, image domain |
| 4 | 2024-11-19 | review | Critical Care | Gaetano Scaramuzzo et al | Electrical impedance tomography monitoring in adult ICU patients: state-of-the-art, recommendations for standardized acquisition, processing, and clinical use, and future directions | 10.1186/s13054-024-05173-x | exclude: ICU EIT consensus review; no decomposition comparison |
| 5 | 2025-01-01 | conference-paper | Lecture notes in electrical engineering | Guobin Gao et al | Variational Mode Decomposition is Applied to the Extraction of Pulmonary Perfusion Signals in Electrical Impedance Tomography | 10.1007/978-981-96-1395-3_58 | exclude: EIT perfusion; variational mode decomposition; one method; conference report of the same work as the Gao et al 2025 journal paper (same group) |
| 6 | 2025-12-16 | article | Measurement Science and Technology | Marlin Ramadhan Baidillah et al | A review of cardiac functional monitoring based on electrical impedance tomography | 10.1088/1361-6501/ae2d86 | exclude: review of cardiac monitoring by EIT; no new comparison |
| 7 | 2026-01-01 | conference-paper | Lecture notes in electrical engineering | Liping Wei et al | Breath-Hold-Free EIT for Pulmonary Perfusion Extraction and Imaging Using Ensemble Empirical Mode Decomposition | 10.1007/978-981-95-6837-6_49 | exclude: EIT perfusion; ensemble EMD without breath-hold; one method (title; no abstract indexed) |
| 8 | 2026-01-09 | review | Physiological Measurement | Yuqiao Peng et al | Electrical impedance tomography for stroke volume monitoring: a narrative review on signal processing, experimental and clinical applications | 10.1088/1361-6579/ae365d | exclude: narrative review, stroke volume by EIT; no new comparison |
| 9 | 2026-04-17 | article | Journal of Nondestructive Evaluation | Laura Homa et al | Multi-Modal Source Separation for Electrical Impedance Tomography | 10.1007/s10921-026-01363-5 | exclude: EIT for aerospace composites (damage vs strain); not physiological |
| 10 | 2026-09-08 **(after 24 Aug 2026)** | article | Biomedical Physics & Engineering Express | Maulina Nurul Nurul Hidayah et al | Cardiac and respiratory signal separation in electrical impedance tomography (EIT) | 10.1088/2057-1976/aea425 | exclude: EIT, simulated thorax: four separation methods (CWT, EMD, ICA, PCA) compared across measurement patterns, scored by correlation of extracted time signals with reference signals and NRMSE of the reconstructed cardiac waveform (abstract). Published 2026-09-08, after the 24 August 2026 search. From the abstract it separates time signals across frames and does not form a single-frame additive contract or the image of the removed component; the full text was not read. Recorded here because it compares several methods in EIT; it does not change the dated statement |

## 3. What this record supports

As of 24 August 2026, in PubMed and in the OpenAlex forward citations of references 19 and 21, DAGAF had not been applied to EIT and no two single-frame measurement-domain decompositions had been compared under one contract. The re-run of 31 October 2026 returns one new work in EIT (Hidayah et al 2026, *Biomed Phys Eng Express*), which compares four separation methods on simulated time signals across measurement patterns; from its abstract it is not a single-frame measurement-domain comparison under one contract, and it postdates the search. The manuscript's statement is dated and unchanged. The limitation stated in the manuscript stands: PubMed does not index the journals in which the DAGAF papers appear, so the OpenAlex citation screen carries most of the evidence, and Scopus and Web of Science were not searched.

Screening was done on titles, venues and, where needed, abstracts (OpenAlex and PubMed); no full text was read for this record. Decisions were made by the authors with AI assistance (Claude, Anthropic), as disclosed in the manuscript's Acknowledgments.
