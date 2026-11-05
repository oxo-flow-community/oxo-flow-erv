oxo-flow-erv
Copyright (c) 2026 oxo-flow-community

This pipeline is a port of SunPast/ERV-pipeline
(https://github.com/SunPast/ERV-pipeline), commit
02ae7481ff8e967ee6532830d0060547b967e57d (main), authored by SunPast.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

---------------------------------------------------------------------
Upstream license

This port is derived from SunPast/ERV-pipeline under the MIT license. The
upstream LICENSE is included verbatim in this repository at
LICENSE.upstream.

The upstream project aligned its FASTQ stage with
OncoHarmony-Network/TCASIA_pipeline (https://github.com/OncoHarmony-Network/TCASIA_pipeline,
MIT) and takes its ERV annotation from the UCSC RepeatMasker table.
No upstream source file is copied verbatim except the license.

---------------------------------------------------------------------
Method citations

This workflow does not introduce a new quantification method; it runs four
published tools on one reference universe. Cite the tool(s) whose output you
use:

- ERVmap — Tokuyama M, Kong Y, Song E, Jayewickreme T, Kang I, Iwasaki A.
  ERVmap analysis reveals genome-wide transcription of human endogenous
  retroviruses. PNAS 2018;115(50):12565-12572.
- Telescope — Bendall ML, de Mulder M, Iniguez LP, et al. Telescope:
  characterization of the transcriptome through annotation of transposable
  elements. Genome Biology 2019;20:213.
- SQuIRE — Yang WR, Ardeljan D, Pacyna CN, Payer LM, Burns KH. SQuIRE reveals
  locus-specific regulation of interspersed repeat expression. Nucleic Acids
  Research 2019;47(14):e79.
- TEtranscripts — Jin Y, Tam OH, Paniagua E, Hammell M. TEtranscripts: a
  package for including transposable elements in differential expression
  analysis of RNA-seq datasets. Bioinformatics 2015;31(22):3593-3599.
