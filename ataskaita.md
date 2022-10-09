# Ataskaita

## Matricos atstumo funkcija

aprašykite, kaip skaičiavote atstumo funkciją;

Atstumo funkcija skirta įvertinti skirtumą tarp dviejų virusų buvo skaičiuojama remiantis kodonų (ir atitinkamai dikodonų) dažnių lentele. Atstumui skaičiuoti buvo naudojama Euklidinis atstumas (https://en.wikipedia.org/wiki/Euclidean_distance) pasitelkus python `scipy` bibliotekos funkciją `pdist`, kuri paskaičiuoja euklidinį astumą tarp dviejų taškų (šiuo atveju virusų) n-matėje erdvėje (https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html). 

Pritaikius python `scipy` bibliotekos funkciją `squareform`, `pdist` funkcijos rezultatai pavaizduoti atstumų matricos pavidalu.

Gautos dvi atstumų matricos - iš kodonų ir dikodonų dažnių lentelės. 

<!-- kokia butu funkcija kuri apskaiciuoja atstuma tarp "sio" viruso ir bet kurio kito viruso? Kaip ivertint kaip jie skirias vienas kitos atzvilgiais. Moduli, dar kanors atimt, arba pakelt kvadratu istraukt sakni, belekaip svarbu kad butu skirtumas -->

#


## Kodonų ir dikodonų medžiai 

Kodonų ir dikodonų medžiai buvo generuojami naudojant *neighbour joining* metodą (http://www.trex.uqam.ca/index.php?action=trex&menuD=1&method=2) bei remiantis anksčiau apskaičiuotomis atstumų matricomis. 

kokie medžiai gavosi su kodonais ir dikodonais;

![alt text](/assets/codon_freq_tree.jpg "Title")

#


## Kodonų ir dikodonų dažnis tarp žinduolių ir bakterijų virusų

Kodonų ir dikodonų dažnis tarp žinduolių ir bakterijų virusų. 

Žvelgiant į gautus virusų medžius, virusai klasterizuojasi grupėmis (obviously)

ar skiriasi kodonų ir dikodonų dažnis tarp žinduolių ir bakterijų virusų, kaip klasterizuojasi virusai. Gal kažkuris virusas labai išsiskyrė? Kokie kodonai/dikodonai labiausiai varijuoja?




Remiantis kodonų lentele, labiausiai varijuojantys kodonai yra , ir.

Remiantis dikodonų lentele, labiausiai varijuojantys dikodonai yra , ir.