import requests
from bs4 import BeautifulSoup

busco = {}
checkm = {}


checkm['ranks'] = ['life', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

checkm['life'] = ['Prokaryote']

checkm['domain'] = ['Archaea', 'Bacteria']

checkm['phylum'] = ['Acidobacteria', 'Actinobacteria', 'Aquificae', 'Bacteroidetes',
       'Chlamydiae', 'Chlorobi', 'Chloroflexi', 'Crenarchaeota',
       'Cyanobacteria', 'Deferribacteres', 'Deinococcus-Thermus',
       'Dictyoglomi', 'Euryarchaeota', 'Firmicutes', 'Fusobacteria',
       'Ignavibacteriae', 'Nitrospirae', 'Planctomycetes',
       'Proteobacteria', 'Spirochaetes', 'Synergistetes', 'Tenericutes',
       'Thaumarchaeota', 'Thermodesulfobacteria', 'Thermotogae',
       'Verrucomicrobia']

checkm['class'] = ['Acidobacteriia', 'Aciduliprofundum', 'Actinobacteria',
       'Alphaproteobacteria', 'Aquificae', 'Archaeoglobi', 'Bacilli',
       'Bacteroidetes Order II. Incertae sedis', 'Bacteroidia',
       'Betaproteobacteria', 'Chlamydiia', 'Chlorobia', 'Chloroflexi',
       'Chroococcales', 'Clostridia', 'Cytophagia', 'Deferribacteres',
       'Dehalococcoidetes', 'Deinococci', 'Deltaproteobacteria',
       'Dictyoglomia', 'Epsilonproteobacteria', 'Erysipelotrichi',
       'Flavobacteriia', 'Fusobacteriia', 'Gammaproteobacteria',
       'Halobacteria', 'Holophagae', 'Ignavibacteria', 'Methanobacteria',
       'Methanococci', 'Methanomicrobia', 'Mollicutes', 'Negativicutes',
       'Nitrosopumilales', 'Nitrospira', 'Nostocales', 'Opitutae',
       'Oscillatoriales', 'Planctomycetia', 'Pleurocapsales',
       'Prochlorales', 'Solirubrobacterales', 'Sphingobacteriia',
       'Spirochaetia', 'Stigonematales', 'Synergistia', 'Thermococci',
       'Thermodesulfobacteria', 'Thermomicrobia', 'Thermoplasmata',
       'Thermoprotei', 'Thermotogae', 'Verrucomicrobiae',
       'Zetaproteobacteria']

checkm['order'] = ['Acholeplasmatales', 'Acidilobales', 'Acidithiobacillales',
       'Acidobacteriales', 'Aciduliprofundum', 'Actinomycetales',
       'Aeromonadales', 'Alteromonadales', 'Aquificales',
       'Archaeoglobales', 'Bacillales', 'Bacteroidales',
       'Bacteroidetes Order II. Incertae sedis', 'Bankia',
       'Bdellovibrionales', 'Bifidobacteriales', 'Burkholderiales',
       'Campylobacterales', 'Cardiobacteriales', 'Caulobacterales',
       'Chlamydiales', 'Chlorobiales', 'Chloroflexales', 'Chromatiales',
       'Chroococcales', 'Clostridiales', 'Coriobacteriales',
       'Cytophagales', 'Deferribacterales', 'Dehalococcoidales',
       'Deinococcales', 'Desulfobacterales', 'Desulfovibrionales',
       'Desulfurellales', 'Desulfurococcales', 'Desulfuromonadales',
       'Dictyoglomales', 'Enterobacteriales', 'Entomoplasmatales',
       'Erysipelotrichales', 'Flavobacteriales', 'Fusobacteriales',
       'Halanaerobiales', 'Halobacteriales', 'Holophagales',
       'Hydrogenophilales', 'Ignavibacteriales', 'Lactobacillales',
       'Legionellales', 'Mariprofundales', 'Methanobacteriales',
       'Methanocellales', 'Methanococcales', 'Methanomicrobiales',
       'Methanosarcinales', 'Methylococcales', 'Methylophilales',
       'Mycoplasmatales', 'Myxococcales', 'Nautiliales', 'Neisseriales',
       'Nitrosomonadales', 'Nitrosopumilales', 'Nitrospirales',
       'Nostocales', 'Oceanospirillales', 'Opitutales', 'Oscillatoriales',
       'Pasteurellales', 'Planctomycetales', 'Pleurocapsales',
       'Prochlorales', 'Pseudomonadales', 'Rhizobiales',
       'Rhodobacterales', 'Rhodocyclales', 'Rhodospirillales',
       'Rickettsiales', 'Rubrobacterales', 'Selenomonadales',
       'Solirubrobacterales', 'Sphingobacteriales', 'Sphingomonadales',
       'Spirochaetales', 'Stigonematales', 'Sulfolobales',
       'Synergistales', 'Thermales', 'Thermoanaerobacterales',
       'Thermococcales', 'Thermodesulfobacteriales', 'Thermoplasmatales',
       'Thermoproteales', 'Thermotogales', 'Thiotrichales',
       'Verrucomicrobiales', 'Vibrionales', 'Xanthomonadales']

checkm['family'] = ['Acetobacteraceae', 'Acholeplasmataceae', 'Acidaminococcaceae',
       'Acidilobaceae', 'Acidithiobacillaceae', 'Acidobacteriaceae',
       'Aciduliprofundum', 'Actinomycetaceae', 'Actinopolysporaceae',
       'Aerococcaceae', 'Aeromonadaceae', 'Alcaligenaceae',
       'Alcanivoracaceae', 'Alicyclobacillaceae', 'Alteromonadaceae',
       'Anaplasmataceae', 'Aquificaceae', 'Archaeoglobaceae',
       'Arthrospira', 'Aurantimonadaceae', 'Bacillaceae',
       'Bacteroidaceae', 'Bankia', 'Bartonellaceae', 'Bdellovibrionaceae',
       'Beijerinckiaceae', 'Bifidobacteriaceae', 'Blattabacteriaceae',
       'Brachyspiraceae', 'Bradyrhizobiaceae', 'Brevibacteriaceae',
       'Brucellaceae', 'Burkholderiaceae', 'Campylobacteraceae',
       'Cardiobacteriaceae', 'Carnobacteriaceae', 'Caulobacteraceae',
       'Cellulomonadaceae', 'Chitinophagaceae', 'Chlamydiaceae',
       'Chlorobiaceae', 'Chloroflexaceae', 'Chromatiaceae',
       'Chroococcidiopsis', 'Clostridiaceae',
       'Clostridiales Family XI. Incertae Sedis',
       'Clostridiales Family XVII. Incertae Sedis', 'Colwelliaceae',
       'Comamonadaceae', 'Conexibacteraceae', 'Coriobacteriaceae',
       'Corynebacteriaceae', 'Coxiellaceae', 'Cryomorphaceae',
       'Cyanobacterium', 'Cyanobium', 'Cyanothece', 'Cyclobacteriaceae',
       'Cytophagaceae', 'Deferribacteraceae', 'Dehalococcoidaceae',
       'Deinococcaceae', 'Dermabacteraceae', 'Dermacoccaceae',
       'Dermatophilaceae', 'Desulfobacteraceae', 'Desulfobulbaceae',
       'Desulfohalobiaceae', 'Desulfomicrobiaceae', 'Desulfovibrionaceae',
       'Desulfurellaceae', 'Desulfurobacteriaceae', 'Desulfurococcaceae',
       'Dictyoglomaceae', 'Ectothiorhodospiraceae', 'Enterobacteriaceae',
       'Enterococcaceae', 'Entomoplasmataceae', 'Erysipelotrichaceae',
       'Erythrobacteraceae', 'Eubacteriaceae', 'Exiguobacterium',
       'Fangia', 'Ferrimonadaceae', 'Fischerella', 'Flammeovirgaceae',
       'Flavobacteriaceae', 'Francisellaceae', 'Frankiaceae',
       'Fusobacteriaceae', 'Gemella', 'Geobacteraceae',
       'Geodermatophilaceae', 'Gloeocapsa', 'Glycomycetaceae',
       'Gordoniaceae', 'Hahellaceae', 'Halanaerobiaceae',
       'Halobacteriaceae', 'Halobacteroidaceae', 'Halomonadaceae',
       'Helicobacteraceae', 'Holophagaceae', 'Hydrogenophilaceae',
       'Hydrogenothermaceae', 'Hyphomicrobiaceae', 'Hyphomonadaceae',
       'Idiomarinaceae', 'Intrasporangiaceae', 'Jonesiaceae',
       'Lachnospiraceae', 'Lactobacillaceae', 'Legionellaceae',
       'Leptolyngbya', 'Leptospiraceae', 'Leptotrichiaceae',
       'Leuconostocaceae', 'Listeriaceae', 'Lyngbya', 'Mariprofundaceae',
       'Methanobacteriaceae', 'Methanocaldococcaceae', 'Methanocellaceae',
       'Methanococcaceae', 'Methanocorpusculaceae', 'Methanomicrobiaceae',
       'Methanoregulaceae', 'Methanosaetaceae', 'Methanosarcinaceae',
       'Methylobacteriaceae', 'Methylococcaceae', 'Methylocystaceae',
       'Methylophilaceae', 'Microbacteriaceae', 'Microchaetaceae',
       'Micrococcaceae', 'Microcoleus', 'Micromonosporaceae',
       'Moraxellaceae', 'Mycobacteriaceae', 'Mycoplasmataceae',
       'Myxococcaceae', 'Nakamurellaceae', 'Nautiliaceae',
       'Neisseriaceae', 'Nitrosomonadaceae', 'Nitrosopumilaceae',
       'Nitrospiraceae', 'Nocardiaceae', 'Nocardioidaceae',
       'Nocardiopsaceae', 'Nostocaceae', 'Oceanospirillaceae',
       'Opitutaceae', 'Oscillatoria', 'Oxalobacteraceae',
       'Paenibacillaceae', 'Parachlamydiaceae', 'Pasteurellaceae',
       'Patulibacteraceae', 'Pelagibacteraceae', 'Pelobacteraceae',
       'Peptococcaceae', 'Peptostreptococcaceae', 'Phyllobacteriaceae',
       'Piscirickettsiaceae', 'Planctomycetaceae', 'Planktothrix',
       'Planococcaceae', 'Pleurocapsa', 'Porphyromonadaceae',
       'Prevotellaceae', 'Prochlorococcaceae', 'Promicromonosporaceae',
       'Propionibacteriaceae', 'Pseudanabaena', 'Pseudoalteromonadaceae',
       'Pseudomonadaceae', 'Pseudonocardiaceae', 'Psychromonadaceae',
       'Pyrodictiaceae', 'Rhizobiaceae', 'Rhodobacteraceae',
       'Rhodobiaceae', 'Rhodocyclaceae', 'Rhodospirillaceae',
       'Rhodothermaceae', 'Rickettsiaceae', 'Rikenellaceae',
       'Rivulariaceae', 'Rubrivivax', 'Rubrobacteraceae',
       'Ruminococcaceae', 'Saprospiraceae', 'Shewanellaceae',
       'Sinobacteraceae', 'Solirubrobacteraceae', 'Sphingobacteriaceae',
       'Sphingomonadaceae', 'Spirochaetaceae', 'Sporolactobacillaceae',
       'Staphylococcaceae', 'Streptococcaceae', 'Streptomycetaceae',
       'Streptosporangiaceae', 'Succinivibrionaceae', 'Sulfolobaceae',
       'Sutterellaceae', 'Synechococcus', 'Synechocystis',
       'Synergistaceae', 'Syntrophomonadaceae', 'Teredinibacter',
       'Thermaceae', 'Thermoactinomycetaceae', 'Thermoanaerobacteraceae',
       'Thermoanaerobacterales Family III. Incertae Sedis',
       'Thermococcaceae', 'Thermodesulfobacteriaceae',
       'Thermodesulfobiaceae', 'Thermomonosporaceae',
       'Thermoplasmataceae', 'Thermoproteaceae', 'Thermotogaceae',
       'Thiomonas', 'Thiotrichaceae', 'Veillonellaceae',
       'Verrucomicrobiaceae', 'Vibrionaceae', 'Xanthobacteraceae',
       'Xanthomonadaceae']

checkm['genus'] = ['Acetobacter', 'Acholeplasma', 'Achromobacter', 'Acidaminococcus',
       'Acidilobus', 'Acidiphilium', 'Acidithiobacillus',
       'Acidobacterium', 'Acidovorax', 'Aciduliprofundum',
       'Acinetobacter', 'Actinobacillus', 'Actinobaculum', 'Actinomadura',
       'Actinomyces', 'Actinoplanes', 'Actinopolyspora', 'Aequorivita',
       'Aerococcus', 'Aeromonas', 'Aggregatibacter', 'Agrobacterium',
       'Agromyces', 'Ahrensia', 'Alcaligenes', 'Alcanivorax',
       'Algoriphagus', 'Aliagarivorans', 'Alicycliphilus',
       'Alicyclobacillus', 'Alishewanella', 'Alistipes', 'Alkaliphilus',
       'Alloscardovia', 'Alteromonas', 'Aminobacterium', 'Amycolatopsis',
       'Anabaena', 'Anaerococcus', 'Anaeromyxobacter', 'Anaplasma',
       'Anoxybacillus', 'Aquimarina', 'Archaeoglobus', 'Arcobacter',
       'Arenibacter', 'Arsenophonus', 'Arthrobacter', 'Arthrospira',
       'Asticcacaulis', 'Atopobium', 'Aurantimonas', 'Azoarcus',
       'Azorhizobium', 'Azospira', 'Azospirillum', 'Azotobacter',
       'Bacillus', 'Bacteroides', 'Bankia', 'Barnesiella', 'Bartonella',
       'Bdellovibrio', 'Bifidobacterium', 'Blastococcus',
       'Blattabacterium', 'Blautia', 'Bordetella', 'Borrelia',
       'Brachyspira', 'Bradyrhizobium', 'Brevibacillus', 'Brevibacterium',
       'Brevundimonas', 'Brucella', 'Buchnera', 'Burkholderia',
       'Butyricimonas', 'Butyrivibrio', 'Caldicellulosiruptor',
       'Calothrix', 'Campylobacter', 'Candidatus Arthromitus',
       'Candidatus Blochmannia', 'Candidatus Liberibacter',
       'Candidatus Pelagibacter', 'Candidatus Pelagibacter-like',
       'Candidatus Phytoplasma', 'Capnocytophaga', 'Carboxydothermus',
       'Cardiobacterium', 'Caulobacter', 'Cellulomonas', 'Cellulophaga',
       'Chitinilyticum', 'Chlamydia', 'Chlamydophila', 'Chlorobaculum',
       'Chlorobium', 'Chloroflexus', 'Chroococcidiopsis',
       'Chryseobacterium', 'Citreicella', 'Citrobacter', 'Citromicrobium',
       'Clavibacter', 'Clostridium', 'Cohnella', 'Collinsella',
       'Colwellia', 'Comamonas', 'Commensalibacter', 'Conchiformibius',
       'Coprobacillus', 'Coprococcus', 'Coprothermobacter',
       'Corynebacterium', 'Coxiella', 'Cronobacter', 'Cupriavidus',
       'Curtobacterium', 'Curvibacter', 'Cyanobacterium', 'Cyanobium',
       'Cyanothece', 'Cycloclasticus', 'Cytophaga', 'Dechloromonas',
       'Dehalobacter', 'Dehalococcoides', 'Deinococcus', 'Delftia',
       'Desulfatibacillum', 'Desulfobacterium', 'Desulfobulbus',
       'Desulfomicrobium', 'Desulfosporosinus', 'Desulfotignum',
       'Desulfotomaculum', 'Desulfovibrio', 'Desulfurococcus',
       'Dialister', 'Dickeya', 'Dictyoglomus', 'Dorea', 'Dyadobacter',
       'Dysgonomonas', 'Echinicola', 'Ectothiorhodospira', 'Edwardsiella',
       'Eggerthella', 'Ehrlichia', 'Elizabethkingia', 'Ensifer',
       'Enterobacter', 'Enterococcus', 'Enterovibrio', 'Entomoplasma',
       'Eremococcus', 'Erwinia', 'Erysipelothrix', 'Erythrobacter',
       'Escherichia', 'Eubacterium', 'Exiguobacterium', 'Facklamia',
       'Faecalibacterium', 'Fangia', 'Ferrimonas', 'Finegoldia',
       'Fischerella', 'Flavobacterium', 'Flexibacter', 'Francisella',
       'Frankia', 'Fusobacterium', 'Gallibacterium', 'Gardnerella',
       'Gemella', 'Geobacillus', 'Geobacter', 'Gillisia', 'Glaciecola',
       'Gloeocapsa', 'Gluconacetobacter', 'Gluconobacter', 'Glycomyces',
       'Gordonia', 'Gramella', 'Granulicatella', 'Granulicella',
       'Haemophilus', 'Hahella', 'Halalkalicoccus', 'Halanaerobium',
       'Haliea', 'Haloarcula', 'Halobacillus', 'Halobacterium',
       'Halobiforma', 'Halococcus', 'Haloferax', 'Halomicrobium',
       'Halomonas', 'Haloquadratum', 'Halorubrum', 'Haloterrigena',
       'Helicobacter', 'Herbaspirillum', 'Hippea', 'Hirschia',
       'Histophilus', 'Hydrogenobacter', 'Hydrogenobaculum',
       'Hymenobacter', 'Hyphomicrobium', 'Idiomarina', 'Isoptericola',
       'Janibacter', 'Janthinobacterium', 'Jeotgalicoccus', 'Jonesia',
       'Jonquetella', 'Kaistia', 'Kandleria', 'Kangiella', 'Kingella',
       'Klebsiella', 'Labrenzia', 'Lachnobacterium', 'Lachnospira',
       'Lactobacillus', 'Lactococcus', 'Laribacter', 'Leeuwenhoekiella',
       'Legionella', 'Leifsonia', 'Leisingera', 'Leptolyngbya',
       'Leptospira', 'Leptotrichia', 'Leucobacter', 'Leuconostoc',
       'Lewinella', 'Listeria', 'Loktanella', 'Luteimonas', 'Lyngbya',
       'Lysinibacillus', 'Lysobacter', 'Magnetospirillum', 'Mannheimia',
       'Maribacter', 'Marinobacter', 'Marinobacterium', 'Marinomonas',
       'Mariprofundus', 'Maritimibacter', 'Massilia', 'Megamonas',
       'Megasphaera', 'Meiothermus', 'Mesoplasma', 'Mesorhizobium',
       'Metallosphaera', 'Methanobacterium', 'Methanobrevibacter',
       'Methanocaldococcus', 'Methanocella', 'Methanococcus',
       'Methanocorpusculum', 'Methanolobus', 'Methanoplanus',
       'Methanosaeta', 'Methanosarcina', 'Methanothermobacter',
       'Methanothermococcus', 'Methanotorris', 'Methylobacter',
       'Methylobacterium', 'Methylococcus', 'Methylocystis',
       'Methylomicrobium', 'Methylomonas', 'Methylophaga',
       'Methylophilus', 'Methylopila', 'Methylosarcina', 'Methylosinus',
       'Methylotenera', 'Methyloversatilis', 'Methylovorus',
       'Microbacterium', 'Microcoleus', 'Micromonospora', 'Microvirga',
       'Mobiluncus', 'Moorella', 'Moraxella', 'Morganella',
       'Mycobacterium', 'Mycoplasma', 'Myroides', 'Natrialba',
       'Natrinema', 'Natronobacterium', 'Natronorubrum', 'Neisseria',
       'Neorickettsia', 'Nesterenkonia', 'Niabella', 'Nisaea',
       'Nitratireductor', 'Nitratiruptor', 'Nitrobacter', 'Nitrosococcus',
       'Nitrosomonas', 'Nitrosospira', 'Nocardia', 'Nocardioides',
       'Nocardiopsis', 'Nostoc', 'Novosphingobium', 'Oceanicaulis',
       'Oceanicola', 'Oceanimonas', 'Oceanobacillus', 'Oceanospirillum',
       'Ochrobactrum', 'Octadecabacter', 'Odoribacter', 'Oenococcus',
       'Oligella', 'Oligotropha', 'Olleya', 'Olsenella', 'Oribacterium',
       'Orientia', 'Ornithobacterium', 'Oscillatoria', 'Oxalobacter',
       'Paenibacillus', 'Pantoea', 'Parabacteroides', 'Parachlamydia',
       'Paracoccus', 'Paraprevotella', 'Parascardovia', 'Pasteurella',
       'Patulibacter', 'Pectobacterium', 'Pediococcus', 'Pedobacter',
       'Pelobacter', 'Pelodictyon', 'Pelosinus', 'Peptoniphilus',
       'Peptostreptococcus', 'Persephonella', 'Phaeobacter',
       'Photobacterium', 'Photorhabdus', 'Planctomyces', 'Planktothrix',
       'Planococcus', 'Pleomorphomonas', 'Pleurocapsa', 'Polaribacter',
       'Polaromonas', 'Polynucleobacter', 'Pontibacillus', 'Pontibacter',
       'Porphyromonas', 'Prevotella', 'Prochlorococcus',
       'Promicromonospora', 'Propionibacterium', 'Propionimicrobium',
       'Proteus', 'Providencia', 'Pseudanabaena', 'Pseudoalteromonas',
       'Pseudobutyrivibrio', 'Pseudoclavibacter', 'Pseudogulbenkiania',
       'Pseudomonas', 'Pseudonocardia', 'Pseudovibrio',
       'Pseudoxanthomonas', 'Psychrobacter', 'Psychroflexus',
       'Psychromonas', 'Pyrobaculum', 'Pyrococcus', 'Rahnella',
       'Ralstonia', 'Rheinheimera', 'Rhizobium', 'Rhodanobacter',
       'Rhodobacter', 'Rhodococcus', 'Rhodonellum', 'Rhodopirellula',
       'Rhodopseudomonas', 'Rhodospirillum', 'Rhodothermus', 'Rickettsia',
       'Riemerella', 'Roseburia', 'Roseiflexus', 'Roseobacter',
       'Roseomonas', 'Roseovarius', 'Rothia', 'Rubrivivax', 'Rubrobacter',
       'Ruegeria', 'Ruminococcus', 'Runella', 'Saccharibacillus',
       'Saccharomonospora', 'Salinicoccus', 'Salinimicrobium',
       'Salinispora', 'Salmonella', 'Saprospira', 'Sediminibacterium',
       'Selenomonas', 'Serratia', 'Shewanella', 'Shigella', 'Slackia',
       'Solirubrobacter', 'Solobacterium', 'Sphaerochaeta',
       'Sphingobacterium', 'Sphingobium', 'Sphingomonas', 'Sphingopyxis',
       'Spirochaeta', 'Spirosoma', 'Sporosarcina', 'Staphylococcus',
       'Staphylothermus', 'Stenotrophomonas', 'Streptococcus',
       'Streptomyces', 'Sulfitobacter', 'Sulfobacillus', 'Sulfolobus',
       'Sulfurihydrogenibium', 'Sulfurimonas', 'Sulfurospirillum',
       'Sutterella', 'Synechococcus', 'Synechocystis', 'Tannerella',
       'Taylorella', 'Tenacibaculum', 'Teredinibacter', 'Tetragenococcus',
       'Thalassospira', 'Thauera', 'Thermacetogenium', 'Thermaerobacter',
       'Thermanaerovibrio', 'Thermoanaerobacter', 'Thermoanaerobacterium',
       'Thermobifida', 'Thermococcus', 'Thermocrinis', 'Thermocrispum',
       'Thermodesulfatator', 'Thermodesulfobacterium',
       'Thermodesulfovibrio', 'Thermogladius', 'Thermoplasma',
       'Thermoproteus', 'Thermosipho', 'Thermotoga', 'Thermus',
       'Thioalkalimicrobium', 'Thioalkalivibrio', 'Thiobacillus',
       'Thiomicrospira', 'Thiomonas', 'Thiothrix', 'Tolumonas',
       'Treponema', 'Turicibacter', 'Ureaplasma', 'Variovorax',
       'Veillonella', 'Verrucomicrobium', 'Vibrio', 'Vulcanisaeta',
       'Weissella', 'Wigglesworthia', 'Wohlfahrtiimonas', 'Wolbachia',
       'Xanthobacter', 'Xanthomonas', 'Xenorhabdus', 'Xylella',
       'Yersinia', 'Zymomonas']


checkm['species'] = ['Acetobacter pasteurianus',
     'Achromobacter piechaudii',
     'Achromobacter xylosoxidans',
     'Acidaminococcus intestini',
     'Acidithiobacillus caldus',
     'Acidithiobacillus ferrooxidans',
     'Acidovorax avenae',
     'Acinetobacter baumanni',
     'Acinetobacter baumannii',
     'Acinetobacter baylyi',
     'Acinetobacter calcoaceticus',
     'Acinetobacter johnsonii',
     'Acinetobacter junii',
     'Acinetobacter lwoffii',
     'Acinetobacter nosocomialis',
     'Acinetobacter radioresistens',
     'Actinobacillus pleuropneumoniae',
     'Actinobaculum schaalii',
     'Actinomyces graevenitzii',
     'Actinomyces naeslundii',
     'Actinomyces neuii',
     'Actinomyces odontolyticus',
     'Aeromonas hydrophila',
     'Aeromonas veronii',
     'Aggregatibacter actinomycetemcomitans',
    'Agrobacterium tumefaciens',
    'Alicycliphilus denitrificans',
    'Alicyclobacillus acidocaldarius',
    'Alloscardovia omnicolens',
    'Alteromonas macleodii',
    'Amycolatopsis mediterranei',
    'Anabaena circinalis',
    'Anaerococcus prevotii',
    'Anaeromyxobacter dehalogenans',
    'Anoxybacillus flavithermus',
    'Archaeoglobus fulgidus',
    'Arcobacter butzleri',
    'Arthrobacter nicotinovorans',
    'Arthrospira platensis',
    'Atopobium vaginae',
    'Azotobacter vinelandii',
    'Bacillus amyloliquefaciens',
    'Bacillus anthracis',
    'Bacillus cereus',
    'Bacillus coagulans',
    'Bacillus licheniformis',
    'Bacillus megaterium',
    'Bacillus methanolicus',
    'Bacillus pumilus',
    'Bacillus subtilis',
    'Bacillus thuringiensis',
    'Bacteroides caccae',
    'Bacteroides dorei',
    'Bacteroides eggerthii',
    'Bacteroides fragilis',
    'Bacteroides massiliensis',
    'Bacteroides ovatus',
    'Bacteroides stercoris',
    'Bacteroides thetaiotaomicron',
    'Bacteroides uniformis',
    'Bacteroides vulgatus',
    'Bacteroides xylanisolvens',
    'Bankia setacea',
    'Bartonella clarridgeiae',
    'Bartonella doshiae',
    'Bartonella elizabethae',
    'Bartonella grahamii',
    'Bartonella quintana',
    'Bartonella tamiae',
    'Bartonella vinsonii',
    'Bartonella washoensis',
    'Bdellovibrio bacteriovorus',
    'Bifidobacterium adolescentis',
    'Bifidobacterium animalis',
    'Bifidobacterium bifidum',
    'Bifidobacterium breve',
    'Bifidobacterium dentium',
    'Bifidobacterium longum',
    'Bifidobacterium pseudolongum',
    'Blautia producta',
    'Blautia wexlerae',
    'Bordetella holmesii',
    'Bordetella pertussis',
    'Borrelia afzelii',
    'Borrelia burgdorferi',
    'Borrelia garinii',
    'Brachyspira hyodysenteriae',
    'Brachyspira pilosicoli',
    'Bradyrhizobium elkanii',
    'Bradyrhizobium japonicum',
    'Brevibacillus laterosporus',
    'Brucella abortus',
    'Brucella canis',
    'Brucella ceti',
    'Brucella melitensis',
    'Brucella ovis',
    'Brucella pinnipedialis',
    'Brucella suis',
    'Buchnera aphidicola',
    'Burkholderia cenocepacia',
    'Burkholderia cepacia',
    'Burkholderia mallei',
    'Burkholderia mimosarum',
    'Burkholderia multivorans',
    'Burkholderia pseudomallei',
    'Butyrivibrio fibrisolvens',
    'Butyrivibrio proteoclasticus',
    'Caldicellulosiruptor kristjanssonii',
    'Campylobacter coli',
    'Campylobacter curvus',
    'Campylobacter fetus',
    'Campylobacter jejuni',
    'Campylobacter showae',
    'Campylobacter upsaliensis',
    'Candidatus Pelagibacter ubique',
    'Candidatus Pelagibacter-like SAR11',
    'Capnocytophaga ochracea',
    'Caulobacter crescentus',
    'Cellulomonas flavigena',
    'Chlamydia muridarum',
    'Chlamydia psittaci',
    'Chlamydia trachomatis',
    'Chlamydophila pneumoniae',
    'Chlamydophila psittaci',
    'Chlorobium phaeobacteroides',
    'Citrobacter freundii',
    'Clavibacter michiganensis',
    'Clostridium acetobutylicum',
    'Clostridium beijerinckii',
    'Clostridium bolteae',
    'Clostridium botulinum',
    'Clostridium butyricum',
    'Clostridium clariflavum',
    'Clostridium clostridioforme',
    'Clostridium difficile',
    'Clostridium glycolicum',
    'Clostridium kluyveri',
    'Clostridium pasteurianum',
    'Clostridium perfringens',
    'Clostridium sporogenes',
    'Clostridium symbiosum',
    'Clostridium thermocellum',
    'Comamonas testosteroni',
    'Corynebacterium accolens',
    'Corynebacterium aurimucosum',
    'Corynebacterium diphtheriae',
    'Corynebacterium efficiens',
    'Corynebacterium glucuronolyticum',
    'Corynebacterium glutamicum',
    'Corynebacterium halotolerans',
    'Corynebacterium jeikeium',
    'Corynebacterium matruchotii',
    'Corynebacterium pseudotuberculosis',
    'Corynebacterium ulcerans',
    'Corynebacterium urealyticum',
    'Coxiella burnetii',
    'Cupriavidus metallidurans',
    'Cupriavidus necator',
    'Cupriavidus taiwanensis',
    'Dechloromonas agitata',
    'Dehalococcoides mccartyi',
    'Deinococcus radiodurans',
    'Delftia acidovorans',
    'Desulfovibrio africanus',
    'Desulfovibrio alaskensis',
    'Desulfovibrio desulfuricans',
    'Desulfovibrio vulgaris',
    'Dickeya dadantii',
    'Dorea longicatena',
    'Edwardsiella tarda',
    'Ehrlichia ruminantium',
    'Elizabethkingia anophelis',
    'Ensifer fredii',
    'Ensifer medicae',
    'Ensifer meliloti',
    'Enterobacter aerogenes',
    'Enterobacter cloacae',
    'Enterococcus casseliflavus',
    'Enterococcus faecalis',
    'Enterococcus faecium',
    'Enterococcus gallinarum',
    'Eremococcus coleocola',
    'Erwinia amylovora',
    'Erwinia pyrifoliae',
    'Escherichia coli',
    'Escherichia fergusonii',
    'Eubacterium cellulosolvens',
    'Eubacterium saburreum',
    'Exiguobacterium sibiricum',
    'Exiguobacterium undae',
    'Facklamia hominis',
    'Finegoldia magna',
    'Francisella noatunensis',
    'Francisella noatunensis noatunensis',
    'Francisella novicida',
    'Francisella philomiragia',
    'Francisella tularensis',
    'Fusobacterium necrophorum',
    'Fusobacterium nucleatum',
    'Gallibacterium anatis',
    'Gardnerella vaginalis',
    'Gemella haemolysans',
    'Geobacillus debilis',
    'Geobacter sulfurreducens',
    'Glaciecola agarilytica',
    'Gluconacetobacter diazotrophicus',
    'Gordonia amicalis',
    'Gordonia polyisoprenivorans',
    'Haemophilus ducreyi',
    'Haemophilus haemolyticus',
    'Haemophilus influenzae',
    'Haemophilus parainfluenzae',
    'Haemophilus parasuis',
    'Halalkalicoccus jeotgali',
    'Haloarcula californiae',
    'Haloarcula vallismortis',
    'Halobiforma lacisalsi',
    'Haloferax denitrificans',
    'Haloferax mucosum',
    'Haloferax sulfurifontis',
    'Haloferax volcanii',
    'Haloquadratum walsbyi',
    'Helicobacter cetorum',
    'Helicobacter pylori',
    'Histophilus somni',
    'Hydrogenobacter thermophilus',
    'Hyphomicrobium denitrificans',
    'Isoptericola variabilis',
    'Jonquetella anthropi',
    'Kandleria vitulina',
    'Klebsiella oxytoca',
    'Klebsiella pneumoniae',
    'Lachnobacterium bovis',
    'Lachnospira multipara',
    'Lactobacillus acidophilus',
    'Lactobacillus amylovorus',
    'Lactobacillus brevis',
    'Lactobacillus buchneri',
    'Lactobacillus casei',
    'Lactobacillus crispatus',
    'Lactobacillus delbrueckii',
    'Lactobacillus fermentum',
    'Lactobacillus gasseri',
    'Lactobacillus helveticus',
    'Lactobacillus iners',
    'Lactobacillus jensenii',
    'Lactobacillus johnsonii',
    'Lactobacillus paracasei',
    'Lactobacillus reuteri',
    'Lactobacillus rhamnosus',
    'Lactobacillus ruminis',
    'Lactobacillus salivarius',
    'Lactococcus garvieae',
    'Lactococcus lactis',
    'Laribacter hongkongensis',
    'Legionella longbeachae',
    'Legionella pneumophila',
    'Legionella sainthelensi',
    'Leptospira meyeri',
    'Leptotrichia goodfellowii',
    'Leuconostoc citreum',
    'Leuconostoc gelidum',
    'Leuconostoc kimchii',
    'Leuconostoc mesenteroides',
    'Listeria monocytogenes',
    'Loktanella vestfoldensis',
    'Lysinibacillus fusiformis',
    'Mannheimia haemolytica',
    'Mariprofundus ferrooxydans',
    'Megasphaera elsdenii',
    'Megasphaera genomosp.',
    'Mesorhizobium ciceri',
    'Mesorhizobium loti',
    'Methanobrevibacter smithii',
    'Methanococcus maripaludis',
    'Methylobacterium extorquens',
    'Methylococcus capsulatus',
    'Methylotenera mobilis',
    'Methylotenera versatilis',
    'Methyloversatilis universalis',
    'Mobiluncus curtisii',
    'Mobiluncus mulieris',
    'Moraxella catarrhalis',
    'Mycobacterium abscessus',
    'Mycobacterium avium',
    'Mycobacterium bovis',
    'Mycobacterium canettii',
    'Mycobacterium gilvum',
    'Mycobacterium intracellulare',
    'Mycobacterium leprae',
    'Mycobacterium marinum',
    'Mycobacterium massiliense',
    'Mycobacterium rhodesiae',
    'Mycobacterium smegmatis',
    'Mycobacterium tuberculosis',
    'Mycoplasma agalactiae',
    'Mycoplasma arginini',
    'Mycoplasma bovis',
    'Mycoplasma canis',
    'Mycoplasma fermentans',
    'Mycoplasma gallisepticum',
    'Mycoplasma genitalium',
    'Mycoplasma hominis',
    'Mycoplasma hyopneumoniae',
    'Mycoplasma hyorhinis',
    'Mycoplasma leachii',
    'Mycoplasma mycoides',
    'Mycoplasma pneumoniae',
    'Mycoplasma synoviae',
    'Myroides odoratimimus',
    'Natronobacterium gregoryi',
    'Natronorubrum tibetense',
    'Neisseria flavescens',
    'Neisseria gonorrhoeae',
    'Neisseria meningitidis',
    'Novosphingobium nitrogenifigens',
    'Oceanicaulis alexandrii',
    'Ochrobactrum intermedium',
    'Oenococcus oeni',
    'Oligotropha carboxidovorans',
    'Orientia tsutsugamushi',
    'Ornithobacterium rhinotracheale',
    'Oxalobacter formigenes',
    'Paenibacillus alvei',
    'Paenibacillus mucilaginosus',
    'Pantoea ananatis',
    'Parabacteroides merdae',
    'Parachlamydia acanthamoebae',
    'Parascardovia denticolens',
    'Pasteurella multocida',
    'Pectobacterium carotovorum',
    'Pediococcus acidilactici',
    'Pediococcus pentosaceus',
    'Pelosinus fermentans',
    'Peptoniphilus lacrimalis',
    'Peptostreptococcus anaerobius',
    'Persephonella lauensis',
    'Phaeobacter gallaeciensis',
    'Photobacterium profundum',
    'Polynucleobacter necessarius',
    'Porphyromonas asaccharolytica',
    'Porphyromonas gingivalis',
    'Prevotella amnii',
    'Prevotella bivia',
    'Prevotella brevis',
    'Prevotella buccae',
    'Prevotella dentalis',
    'Prevotella denticola',
    'Prevotella intermedia',
    'Prevotella maculosa',
    'Prevotella melaninogenica',
    'Prevotella micans',
    'Prevotella oralis',
    'Prevotella oris',
    'Prevotella ruminicola',
    'Prevotella timonensis',
    'Prevotella veroralis',
    'Prochlorococcus marinus',
    'Propionibacterium acidifaciens',
    'Propionibacterium acidipropionici',
    'Propionibacterium acnes',
    'Propionibacterium avidum',
    'Propionimicrobium lymphophilum',
    'Proteus mirabilis',
    'Providencia alcalifaciens',
    'Providencia stuartii',
    'Pseudoalteromonas haloplanktis',
    'Pseudoalteromonas luteoviolacea',
    'Pseudoalteromonas piscicida',
    'Pseudobutyrivibrio ruminis',
    'Pseudomonas aeruginosa',
    'Pseudomonas brassicacearum',
    'Pseudomonas chlororaphis',
    'Pseudomonas fluorescens',
    'Pseudomonas fulva',
    'Pseudomonas mendocina',
    'Pseudomonas monteilii',
    'Pseudomonas putida',
    'Pseudomonas stutzeri',
    'Pseudomonas syringae',
    'Pseudomonas thermotolerans',
    'Pseudomonas umsongensis',
    'Pseudomonas viridiflava',
    'Pseudoxanthomonas suwonensis',
    'Pyrococcus furiosus',
    'Ralstonia pickettii',
    'Ralstonia solanacearum',
    'Rhizobium etli',
    'Rhizobium leguminosarum',
    'Rhodobacter sphaeroides',
    'Rhodococcus equi',
    'Rhodococcus erythropolis',
    'Rhodococcus pyridinivorans',
    'Rhodonellum psychrophilum',
    'Rhodopseudomonas palustris',
    'Rhodospirillum rubrum',
    'Rhodothermus marinus',
    'Rickettsia canadensis',
    'Rickettsia massiliae',
    'Rickettsia prowazekii',
    'Rickettsia rickettsii',
    'Rickettsia slovaca',
    'Rickettsia typhi',
    'Riemerella anatipestifer',
    'Rothia dentocariosa',
    'Rothia mucilaginosa',
    'Rubrivivax benzoatilyticus',
    'Ruminococcus albus',
    'Ruminococcus flavefaciens',
    'Ruminococcus gnavus',
    'Saccharomonospora azurea',
    'Salinispora arenicola',
    'Salinispora pacifica',
    'Salinispora tropica',
    'Salmonella enterica',
    'Saprospira grandis',
    'Selenomonas artemidis',
    'Selenomonas bovis',
    'Selenomonas ruminantium',
    'Selenomonas sputigena',
    'Serratia marcescens',
    'Serratia plymuthica',
    'Serratia proteamaculans',
    'Shewanella baltica',
    'Shigella boydii',
    'Shigella dysenteriae',
    'Shigella flexneri',
    'Shigella sonnei',
    'Solobacterium moorei',
    'Sphingomonas melonis',
    'Sphingomonas phyllosphaerae',
    'Spirochaeta thermophila',
    'Staphylococcus aureus',
    'Staphylococcus epidermidis',
    'Staphylococcus hominis',
    'Staphylococcus lugdunensis',
    'Staphylococcus pseudintermedius',
    'Stenotrophomonas maltophilia',
    'Streptococcus agalactiae',
    'Streptococcus anginosus',
    'Streptococcus australis',
    'Streptococcus bovis',
    'Streptococcus dysgalactiae',
    'Streptococcus equi',
    'Streptococcus gallolyticus',
    'Streptococcus infantarius',
    'Streptococcus infantis',
    'Streptococcus intermedius',
    'Streptococcus mitis',
    'Streptococcus mutans',
    'Streptococcus oralis',
    'Streptococcus parasanguinis',
    'Streptococcus parauberis',
    'Streptococcus pneumoniae',
    'Streptococcus pseudopneumoniae',
    'Streptococcus pseudoporcinus',
    'Streptococcus pyogenes',
    'Streptococcus salivarius',
    'Streptococcus sanguinis',
    'Streptococcus suis',
    'Streptococcus thermophilus',
    'Streptococcus urinalis',
    'Streptomyces cattleya',
    'Sulfobacillus acidophilus',
    'Sulfolobus acidocaldarius',
    'Sulfolobus islandicus',
    'Sulfolobus solfataricus',
    'Sutterella wadsworthensis',
    'Synechococcus elongatus',
    'Taylorella asinigenitalis',
    'Taylorella equigenitalis',
    'Teredinibacter turnerae',
    'Thauera linaloolentis',
    'Thermacetogenium phaeum',
    'Thermoanaerobacterium thermosaccharolyticum',
    'Thermobifida fusca',
    'Thermus oshimai',
    'Thermus scotoductus',
    'Thermus thermophilus',
    'Thioalkalivibrio thiocyanoxidans',
    'Thiobacillus denitrificans',
    'Treponema denticola',
    'Treponema pallidum',
    'Treponema vincentii',
    'Ureaplasma parvum',
    'Ureaplasma urealyticum',
    'Variovorax paradoxus',
    'Veillonella atypica',
    'Veillonella parvula',
    'Vibrio alginolyticus',
    'Vibrio cholerae',
    'Vibrio fischeri',
    'Vibrio harveyi',
    'Vibrio mimicus',
    'Vibrio parahaemolyticus',
    'Vibrio splendidus',
    'Vibrio vulnificus',
    'Wigglesworthia glossinidia',
    'Wohlfahrtiimonas chitiniclastica',
    'Wolbachia endosymbiont of Culex quinquefasciatus',
    'Xanthomonas axonopodis',
    'Xanthomonas campestris',
    'Xanthomonas oryzae',
    'Xylella fastidiosa',
    'Yersinia enterocolitica',
    'Yersinia pestis',
    'Yersinia pseudotuberculosis',
    'Zymomonas mobilis']


def get_busco_lineages_and_urls(base_url='https://busco-data.ezlab.org/v5/data/lineages/'):
    """
    Fetches and prints the list of files and their URLs from the given URL.

    :param base_url: URL of the page to scrape
    :return: List of tuples containing file names and their corresponding URLs
    """
    response = requests.get(base_url)
    response.raise_for_status()  # Raise an exception for HTTP errors

    # Parse the HTML content with BeautifulSoup
    soup = BeautifulSoup(response.content, 'html.parser')

    # Find all anchor tags (<a>) that have an href attribute
    links = soup.find_all('a', href=True)

    # Extract file names and their full URLs
    files = [(link.get_text().split("_odb10")[0], base_url + link['href']) for link in links if link['href'].
endswith(('tar.gz', 'zip', 'tar', 'gz'))]

    return dict(files)

busco = get_busco_lineages_and_urls()


