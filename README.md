## NEWPACK

Experimental Lossless Compressor

Copyright (c) 2020 David Bryant.

All Rights Reserved.

Distributed under the [BSD Software License](https://github.com/dbry/newpack/blob/master/license.txt).

## What is this?

This is an implementation of an experimental general-purpose lossless data compressor, the genesis of which was the "fast" DSD compression mode I developed for [WavPack](https://github.com/dbry/WavPack). I discovered early on that several regular lossless compressors (e.g. `bzip2`) do a surprisingly decent job on DSD (1-bit PCM) audio files, and after I developed the "fast" DSD mode of WavPack I discovered that it could often do a decent job on other types of data (like text).

This implementation pre-scans each buffer of data to be compressed and creates a probability model for each byte which is based on some number of previous bytes in the stream. This model is compressed with simple RLE schemes and sent to the output file, and then the actual data is encoded using the model with a range coder and appended to the stream.

The model can incorporate from 0 to 8 previous bytes in the prediction, however memory constraints limit the number that are actually practical. For text files around 2-3 previous bytes are generally optimum, and for binary files 1-2 bytes work best. Only if the file uses a very small symbol set can 4 or more bytes be used successfully in the model. Having zero history bytes simply means that no history is considered and only the frequency of each isolated symbol is stored.

## How does it compare?

Because the model can become very large, better compression is achieved with large blocks of data (assuming the data is somewhat homogeneous) and also works better with text files, etc., with smaller symbol sets. It often performs very poorly compared to other compressors, but in some cases compresses better than all compressors I've compared it with. For example, DSD audio files (which makes sense) but also the `Large Canterbury` corpus `E.coli` file:

Compressor | Compressed Size
-----------|----------------
None       | 4.638,690
gzip       | 1,299,066 
bzip2      | 1,251,004
compress   | 1,218,349
[lzw-abi](https://github.com/dbry/lzw-ab)     | 1,215,447 
xz         | 1,186,180
newpack    | 1,126,793

Another file that it does very well on (only losing to bzip2) is the `x-ray` file from the `Silesia` corpus:

Compressor | Compressed Size
-----------|----------------
None       | 8,474,240
compress   | 7,008,473
[lzw-abi](https://github.com/dbry/lzw-ab)     | 6,736,569
gzip       | 6,037,713
xz         | 4,491,264
newpack    | 4,369,618
bzip2      | 4,051,112

## Building

To build the demo app on Linux or OS-X:

> $ gcc -O3 *.c -lm -o newpack

The "help" display from the demo app:

```
 NEWPACK  Experimental Lossless Compressor  Version 0.0.1
 Copyright (c) 2020 David Bryant. All Rights Reserved.

 Usage:   NEWPACK [-options] infile outfile
           specify '-' for stdin or stdout

 Options: -d     = decompress
          -h<n>  = history bytes (compress only, else find best)
          -b<n>  = block size (compress only, default = 16 MB)
          -m     = decode random output based solely on model
                   (works for the first block of file only)
          -vv    = very verbose (include internal details)
          -v     = verbose

 Warning: EXPERIMENTAL - DON'T EVEN THINK OF ARCHIVING WITH THIS!!
```
## Using the model to generate random output

The decompressor includes an option to simply use the model to generate random output (ignoring the actual encoded data) and this can provide an intuitive feel for how the compressor works. For example, I compressed a file containing simply lower-case text from many academic papers utilizing from 0 to 4 history bytes and then generated output from each using only the model. It is obvious how increasing the length of the history produces a better approximation of real text with more actual correctly spelled words. Interestingly, the best compression was achieved with the 4-byte history, even after accounting for the much larger model. These were the results including the size of the stored model:
 ### 0 history bytes, model is 32 bytes in length
> cenlotSrnthbl  esnoodo h emenewtfst esft  osi aconl agnton fetcbo eicieua
> cabnnnllnnm gcoietehharwn pdihntdeesn ifol e estn osnr ronit oet nt thenaac
> ansnyothdesrinaotto ticiettesg adituntvhd Saituun dihur ea tuworramlnt r
> orlawaotndggvwg  tge on p etrii ntdtptauledcbnnr  teo e lmlc,rtdfntp aihutehhec
> taes acin ywe  ewtrea gra  eieee idaantio eorntes ncoyetiirichnfaphsncaettaeoniswcde
> cande,yt rltes solaSosetlcrp esteaaeu anettaboa mtgunati e rn ncr eoihswsrtag
> iee redhfarmoptvgamsrus ta,reo ilne lt  tsubadts uaoea mtnnenws me nar ledp
> tetnslarhmtannillo t i rhecy seh rt tefin eigcmrynnsgsstbhdnarwnseuehsoih
> gcire eocat  dh  ntnatpsodbf aierm c nhooit mv meeesec uioeorsefntie ysvstnyh
> aemeibihsr s ni diarkuma eseiianitnage ensfo rgtt  ii nirt eeiioseereomitgoneloi
### 1 history byte, model is 819 bytes in length
> me is acthec f tist t ivans re whatid watinante tstesthisepe t igatoremeresyphe
> ond atessun aligs d ag s derefemo red habe a gk mpran ay tis ons reathedlo
> fre re po atid aty emanibel tre inendenthofombjute des imuganstqucin mathes
> mpveldbomerthy thexpemed agn w cheaso pezrfor areflysess kive vergre whe a
> lem iabdin t amorsend ofupede findithy witas ing tic almply shalless tion
> onont julenvedengl ansuszatevarel theseroumerassthan acolical latluthelo mowhe
> ofentime sus tists th anticowither d pinisopstornghr le oreder theguans topiofideroun
> fwhalathilens ff in bofordeer ngex s atsisppe aculedg we ft in steifcos thderespe
> tin ase concthedear halp the orontescorofly araticticiore biqubovene cesnives
> ks itolorct le maldiefftharofathecenonalts gope ctuspomabed t chathonthali
### 2 history bytes, model is 13,066 bytes in length
> dapprotival tondent ithavel asizationt in repts then hista of in the beejech
> ing thower prowe dinges offich complideliker ad the how atis des ing the and
> witeratorms a la is ch avirs aftes in ase elat to al of tionauggowellegions
> fory evernial so tjc conatuattemist congsums tic pling tharinalit ard but
> lat allowledis elograch sitiall s thow ary smaliogichis foreglizature anconstal
> rels orove nuffecomphave afecrulat ozon masse pergy ed annothato may a gions to
> thworand ved cormare cousivint pent hmand an clue mose dientation red to ged
> lown forld tonal hericsfors prindificition ounduate res wer vain th s iii
> maignits forructo fj inuffere a worallydatortiest dentempecomy any amal
> manationfach thell goves outione ster ball stral bet tradurights so prom of
### 3 history bytes, model is 84,876 bytes in length
> diffected tood oth social lis acquirectiote inly line and of set it metrich
> detent experal difficially foresou fowledge efferime in betwo aimselem mean
> althostution responent the devery ided in the vehic pation the industoric
> respectices anal and could manages thed numerthe thanguant the changle of
> and their fall appeach or in est regulatice notion not mq sets of be enessess
> the diressagemers and throutions of the firs first thus bund part threemendersburgues
> tituteder uniquence of on both limiledge ther vial in same inbrows therject
> inditory between homke and youndama if wher contate flue infor show homent
> and with fundefinity how these wated the pd freques the proach futual data
> s detee as but the progency molengage bankinduced crost the same bk influidealth
### 4 history bytes, model is 333,545 bytes in length
> different that was ostered with morphic effective successessing a neight of
> the mode while partical to error children on the diversity i would a forge
> potention we housing between showing to ense that in federable frozen fox
> and been in the cannot relate qmax of a seen the lbadha in tm time presearch
> splittle here interacteriod had be tahiti ther not a variable modity to componentity
> tourse case of charact to be shown they s going has a were iliar was it as
> over or industration relation of a posite assy card s ideologisting us there
> thetical approach story were inson where frequence stative democrationalisting
> lemma and regulated rhymic different the sets of added future a performan
> other only it is quantum bickey muslim a generation agree managemandard for

## Future improvements?

There are three areas where further experimentation and improvement are
possible:

1. The models (both the used history bin mask and the probability tables) are compressed with simple RLE algorithms (one based on bits and one on bytes). Using the range coder for these would improve the compression.

2. The way that probabilities are stored is not optimized for files that contain a few very common symbols and many relatively rare symbols (because the dynamic range of the probabilities is only 1:160). This is not very common, but fixing this could greatly improve compression on some files.

3. The compressor can work poorly on heterogeneous files like some tar files or executables. One solution would be to prescan large buffers of data and divide the data based on the symbols used and compress them separately, or even combine similar areas.


