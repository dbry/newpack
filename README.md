## NEWPACK

Experimental Lossless Compressor

Copyright (c) 2020 David Bryant.

All Rights Reserved.

Distributed under the [BSD Software License](https://github.com/dbry/newpack/blob/master/license.txt).

## What is this?

This is an experimental general-purpose lossless data compressor, the genesis of which was the "fast" DSD compression mode developed for [WavPack](https://github.com/dbry/WavPack). I discovered early on that several regular compressors (e.g., `bzip2`) do a surprisingly decent job on DSD (1-bit PCM) audio files. Conversely, after I developed the "fast" DSD mode of WavPack, I discovered that it could often do a decent job on other types of data (like text).

This implementation scans each buffer of data to be compressed and creates a probability model for each byte which is based on the hash of some number of previous bytes in the stream. This model consists of a bitmask representing which hashes are actually present in the file and a symbol probability table for each hash that appeared. These two model components are then compressed recursively and sent to the output file, and then the actual data is encoded using the model with a range coder and appended to the end of the encoded block.

The history depth is controlled with a single parameter that goes from 0 to 7, with zero representing no history (i.e., only the frequency of eash isolated symbol is stored) and 7 representing the maximum practical history size (up to a 30-bit hash and up to 2&#94;20 probability tables). The default is 5 (24-bit hash and 2&#94;16 tables) for the default block size of 20,000,000 bytes.

For version 0.0.2 there are options to preprocess each block for better compression. First, if the data is interleaved or periodic for some other reason, then deinterleaving the data with a given "stride" can often improve performance. For example, 16-bit stereo audio data works much better deinterleaved with a stride of 4. A particular stride may be explicitly specfied with `-i<n>`, or there is an option `-i` that will scan each block in an attempt to detect any periodicity and will check if that produces better compression.

The other preprocessing operation is for numerical data (as opposed to symbolic) and performs an arithmetic "delta" operation for each byte. Again, there is an option to force this operation `-nn` and an option to try the option for each block and use it if it offers improvement `-n`. 

Finally, there is an command-line option`-a` to enable tests for both preprocessors. This is rather slow, but can only improve the compression because if the preprocessors don't actually generate an improvement, they are not used.

**It's important to note that this is *very* experimental. It may not faithfully restore some files (and do so silently). Executables built on different platforms may not be compatible with each other, and future versions will probably not decompress files created with this version. It's probably full of undefined behavior (UB) and potentially exploitable memory corruption bugs and does virtually no sanity checking of incoming data. Certainly not recommended as a day-to-day compressor.**

## How does it compare?

Because the model can become very large, better compression is achieved with large blocks of data (assuming the data is somewhat homogeneous). It often performs very poorly compared to other compressors, but in some cases compresses better than all compressors I've compared it with. For example, DSD audio files (which makes sense, considering its origin) but also other non-text files like the `Large Canterbury` corpus `E.coli` file:


Compressor | Compressed Size
-----------|----------------
None       | 4,638,690
gzip       | 1,299,066 
bzip2      | 1,251,004
compress   | 1,218,349
[lzw-ab](https://github.com/dbry/lzw-ab)     | 1,215,447
xz         | 1,186,180
ppmd       | 1,138,541
newpack    | 1,124,601

Another file it performs the best on is the first billion digits of pi text file found on the web. This file is interesting for two reasons. First, because the digits of pi are uncorrelated we shouldn't benefit from any history and should simply utilize symbol frequency (and we do). Also, because of this, it's possible to accurately calculate the maximum compression possible using the formula `10^9 * log (10) / log (256)` (ignoring the case where the decompressor actually calculates pi):

Compressor | Compressed Size
-----------|----------------
None       | 1,000,000,002
gzip       | 469,005,528
[lzw-ab](https://github.com/dbry/lzw-ab)     | 455,135,941
compress   | 452,298,851
xz         | 436,999,340
bzip2      | 430,858,140
ppmd       | 418,108,436
newpack    | 415,244,883
idealized  | 415,241,012

Newpack is only 3871 bytes longer than the absolute minimum possible, an inflation of less that 1 part in 100,000!

Another file that it betters other compressors (with the -a option which enables both 2-byte deinterleaving and numeric mode) is the `x-ray` file from the `Silesia` corpus:

Compressor | Compressed Size
-----------|----------------
None       | 8,474,240
compress   | 7,008,473
[lzw-ab](https://github.com/dbry/lzw-ab)     | 6,736,569
gzip       | 6,037,713
xz         | 4,491,264
bzip2      | 4,051,112
ppmd       | 4,024,809
newpack    | 3,967,934

Uncompressed images often compress very well with newpack. Here's a personal RGB image `Mendocino.bmp` preprocessed using forced numeric mode (-nn):

Compressor | Compressed Size
-----------|----------------
None       | 6,220,922
compress   | 2,897,271
[lzw-ab](https://github.com/dbry/lzw-ab)     | 2,615,720
gzip       | 2,480,082
xz         | 1,846,188
bzip2      | 1,801,375
ppmd       | 1,768,928
newpack    | 1,555,358

The file `dickens`from the `Silensia` corpus demonstrates the text file performance of newpack (this time depth -6). While this is not shabby (losing only to ppmd), the fact that newpack compresses the entire 10M+ file in a single block gives it a sizable advantage over the others:

Compressor | Compressed Size
-----------|----------------
None       | 10,192,446
compress   | 4,007,235
gzip       | 3,868,778
[lzw-ab](https://github.com/dbry/lzw-ab)     | 3,733,338
xz         | 2,831,212
bzip2      | 2,799,520
newpack    | 2,691,457
ppmd       | 2,497,415

Finally, here's an example file that we do very poorly on (and I'm not really sure why). It's the `xml` file, again from the `Silensia` corpus, and this is compressed with 1M blocks:

Compressor | Compressed Size
-----------|----------------
None       | 5,345,280
compress   | 1,133,617
[lzw-ab](https://github.com/dbry/lzw-ab)     | 1,049,120
newpack    | 785,732
gzip       | 691,992
ppmd       | 604,968
bzip2      | 441,186
xz         | 434,892

## Building

To build the demo app on Linux or OS-X:

> $ gcc -O3 *.c -lm -o newpack

The "help" display from the demo app:

```
 NEWPACK  Experimental General-Purpose Lossless Compressor  Version 0.0.2
 Copyright (c) 2020 David Bryant. All Rights Reserved.

 Usage:   NEWPACK [-d] [-options] infile outfile
           specify '-' for stdin or stdout

 Options: -d     = decompression (default is compression)
          -a     = use all specialty modes (equivalent to -irn)
          -[1-7] = probability model depth (default = 5 for 20 MB block)
          -0     = no history employed in model (symbol frequency only)
          -e     = exhaustive search (very slow for very little gain)
          -i     = automatically detect interleave and use if better
          -i<n>  = force specified interleave stride (1-16, default = 1)
          -b<n>  = specify block size (1 MB - 100 MB, default = 20 MB)
          -l     = use long blocks (100 MB; history depth set to -6)
          -s     = use short blocks (4 MB; history depth set to -4)
          -t     = use tiny blocks (1 MB; history depth set to -3)
          -n     = try numerical data type preprocessing and use if better
          -nn    = force numerical data type preprocessing (i.e., deltas)
          -r     = try simple run-length encoding and use if better
          -rr    = force simple run-length encoding (really just for testing)
          -m     = decode random output based solely on model
          -vv    = very verbose messaging (include internal details)
          -v     = verbose messaging

 Warning: EXPERIMENTAL - DON'T EVEN THINK OF ARCHIVING WITH THIS!!
```
## Using the model to generate random output

The decompressor includes an option to simply use the model to generate random output (ignoring the actual encoded data) and this can provide an intuitive feel for how the compressor works. For example, I compressed a file containing simply lower-case text from many pieces of fiction literature at depths -0, -1, -3, and -5 and then generated output from each using only the model. It is obvious how increasing the history depth produces a better approximation of real text with more actual correctly spelled words. The best compression was achieved with the greatest depth, even after accounting for the much larger model, and here the compressor beats everything I tried against it except ppmd. These were the results including the sizes of the stored modeli and the encoded data:
### depth -0, model = 49 bytes, encoding = 2652832 bytes, total = 2652901 bytes
> mrrmyotnanh  atohoi ef eso nwrr w oeh a ri iiwtnhaetsr  la girner esfa mn
> sp  e odn  soe anedk nt cr oea  r esrnarg o  nl ahaisagtf oeoe rnieoaf rei
> fen  ctso to  uaote t ag   iaekk  ytegr  ofewswvwo a  tlnoniae el aostdehttmatrif
> irreiagelooo  a rfnenp rnle  snialta an sd yveadsrwtoroateta foooentobh srtaht
> s  tlnnnhtooshdn trxbrpmt nawnhotcaacishh io ddwhoefuoihisn nhfnhed poailcniaieaiic
> blrtteoy c ci clsuethoyshewwttheilu hn l nta eaarc idy uh niodne  earnnotdiiny
> eeoi  yemh gu rteaseitnotb heou ddnk sl hhso lwd aee tlk egbrmtishetts  otl
> eh  tinsfh se  du rfhi   ueitdon i heincgrp  b ninrmu eui sacisupsp ednedd
> neh ha irtlty  ukit yisu fmaiseeh nlo ypngwge ntaes eeespe peafha hthiitrfo
> inot  o hhhn v oenrr  toti h tprnnlwaahshee oon o tnioe seihsseta oi pd aceagoefit
### depth -1, model = 5433 bytes, encoding = 1866409 bytes, total = 1871842 bytes
> wrayilee for sorembed ancer hes hater tho mores ther hallat womessing fum
> th as throbdrept his a paill to sigist he book no sm bl ance clecid thalle
> le s an de ang hinditere peclownst yon agrah iuitect maskind evepadulact ing
> thimsect sal the kes eass mat foad tat des thavearthisioned had onsing a vic
> heatent daypterelf moulde he frothath one cick to befess cace bl bacecf at
> se res ting fho and bact thime gis cong bing it boypchis hosto anoughim bon
> ambece fome morm home aftent ited quecamered dol rest eastiongs boncer witud
> mit appp hem awrighent yatimps uorackmout hin fou comm ounich tooks a gor
> suculint of you ther the anyarts ce lrome ted s thand ing himousetalc red
> that pr pen onps ppbd toun my attes al the mellabled tot bruland itede le
### depth -3, model = 31468 bytes, encoding = 1559295 bytes, total = 1590763 bytes
> whate i ree froomming her cruetted liveniting ther have plumines a dempty
> ands done to nt and chall on he augh the surt you girls lity any roying tain
> compace appender was notify sprat ming of no mothe shairl said in his pution
> surpost and fout leat s windonsied wife word i he telists uspent on to if
> you thes the kined goine he ours the merat thoutchild left my dest the sto
> member is stalkints fered aways ing rom me motiona emmark hour aways and hear
> pong miled mirach havering you want obed told are the anybod you wer mot for
> hat peavoice weened shed wilight ted det as is sour what their corrous it
> ins exclim and said down as one vulstears you nt kinge hous bes winkin ent
> forward pute ther mayark anot the but till why fluxurings the back sepighters
### depth -5, model = 154573 bytes, encoding = 1247796 bytes, total = 1402369 bytes
> was it wakan the one i love hung up much in the fade this it i head petter
> that stuck door least dream burth the woman answer the germanshing a children
> patrictly pilot sadness killestic jerking the cans touchest you oh abour poems
> that i certain the erudimense orde a frank buzzing they hair widered to trange
> to would had has beyond starlette fore to give went to speech the knew experiod
> or the breason somethods she identing incapable i wilder a nice to your pall
> sleep from the hallways some tongue holy abrupt this wayward of new bootlegged
> in ther matter was but him to giving agains hike the crowned the ruck speaking
> to tl off his his pudding it an as cave had again in you wanted their little
> food foam out them the tipped to him said she curtly motion of the is a first

For that file, here are the results for other compressors:

Compressor | Compressed Size
-----------|----------------
None       | 5,196,229
gzip       | 1,940,202
[lzw-ab](https://github.com/dbry/lzw-ab)     | 1,877,627
xz         | 1,507,212
bzip2      | 1,470,259
newpack    | 1,402,369
ppmd       | 1,326,482

## Future improvements?

There are at least three areas where further experimentation and improvement are possible:

1. The compression operation is currently very slow because of the extensive searching for the optimum hash configuration and history depth for each block. I believe that this could be improved significantly, and does affect the decompression speed, but it's problematic nevertheless.

2. A significant portion of history hashes often result in only a single output code, but even for these cases we still store an entire probability table (albeit compressed). Optimizing for this case might be more efficient, and would probably result in faster decode.

3. The compressor can work poorly on heterogeneous files like some tar files or executables. One solution would be to prescan large buffers of data and divide the data based on the symbols used and compress them separately, or even combine similar areas.


