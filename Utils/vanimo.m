function cmap = vanimo(N)
% Colour map by Fabio Crameri: 
%   Crameri, F. (2018), Scientific colourmaps, Zenodo, doi:10.5281/zenodo.1243862
%

cmap = ...
  [1 0.80345754277515158 0.992153456360599;
   0.99396554896364175 0.79197062790022321 0.98373971417282413;
   0.98791049723394353 0.78051745089175628 0.9753498657738876;
   0.9818511346768195 0.76910343196256326 0.96699397191514291;
   0.975778198216012 0.75774218138091642 0.95866754777432084;
   0.96970706799259054 0.74642730432947124 0.95037180862849924;
   0.9636311931786794 0.73517284067845456 0.942114358685549;
   0.95755446568097857 0.72397206606607523 0.93388652871246847;
   0.95147333474050177 0.71283563984359133 0.92570475435406252;
   0.94539274642622517 0.70176721328741543 0.91756029264625594;
   0.93930755896397589 0.69076795098703458 0.909445967190424;
   0.93321793051159385 0.67983984379363349 0.90137390391266781;
   0.92713118642585668 0.668991164329494 0.89333797641671431;
   0.92104166776414587 0.65821140221073338 0.88534362617825579;
   0.91495181652434909 0.6475066554521135 0.8773821562250701;
   0.90885708947593968 0.6368878235549259 0.8694611946039692;
   0.902761550193271 0.626340568473632 0.86157572303923091;
   0.89665896280838109 0.615881799366969 0.85372212654883173;
   0.89054769365667819 0.60550501995079886 0.84590949202674093;
   0.88443582479841243 0.59522174175429909 0.8381339183367772;
   0.87831332183115429 0.58502958728762389 0.83038982068471312;
   0.87218638462001563 0.5749149050147252 0.82267923218296757;
   0.86604599536626492 0.56490047171278635 0.81499868528969432;
   0.85989852477089246 0.55498695949964372 0.8073559213963255;
   0.853732255527154 0.54516849715508353 0.79974497045009707;
   0.84756045840022232 0.53544114687670474 0.79215792464704082;
   0.84137685425848618 0.5258278974142635 0.78460986888591133;
   0.83517379037149453 0.51630446776779315 0.77708808520821437;
   0.82895800980628886 0.50689883949344627 0.76958769869497223;
   0.82271950153502893 0.49760808587484429 0.762122991178999;
   0.81646522490379247 0.48841491098802281 0.75469340349070413;
   0.81018422934759327 0.47933965936361062 0.74728263019269292;
   0.80389107637367174 0.47038263926943025 0.739903745893086;
   0.79756948654038817 0.46154321738374582 0.7325477379312535;
   0.79122859232347675 0.45283102049993063 0.7252205833019929;
   0.78486662274805963 0.44424214581369148 0.71791625538751924;
   0.77847299942364523 0.43577863795270649 0.71064093514389881;
   0.77205748404201524 0.42745018123637235 0.703391685590954;
   0.76561824006466861 0.41924907101602482 0.696169746080941;
   0.7591434096128632 0.41118484010596806 0.68897406710778164;
   0.75264416194455286 0.40326677116629878 0.68179136495386539;
   0.74612022655627142 0.39548568510522436 0.67464629144270361;
   0.739569365761163 0.38783436459744841 0.66752261215976794;
   0.73297230930204327 0.38033972509559744 0.66040940284752891;
   0.72633960261321773 0.37296942565116387 0.65331404080912736;
   0.71966618947402639 0.36574698218586166 0.64623316881236914;
   0.71293397380308265 0.35863870467883319 0.63915254519525766;
   0.70615340339595145 0.35166430204081967 0.63206392134160838;
   0.69929003758251151 0.34480657250814983 0.62495542570898477;
   0.69235886994097662 0.33804169044374616 0.61781509267269519;
   0.68531931347438313 0.33137316291102387 0.61063740637985664;
   0.67817391945897687 0.32479409350437438 0.60340060065959844;
   0.67090956469857588 0.31829649787490372 0.59609086312612269;
   0.66351484652468984 0.31184286071176143 0.58869968759093627;
   0.65597573092974526 0.305490221663851 0.581229960052111;
   0.648281238342906 0.29916667114887668 0.57366307854088083;
   0.64044649702792888 0.29288742870407691 0.56599215283911375;
   0.63244939409104406 0.28666960544644987 0.55821595047691;
   0.62429770805507911 0.28051162179889155 0.55034726613316287;
   0.615982490635849 0.27442150609528188 0.54236825094289265;
   0.60752121682249927 0.2683841723199869 0.53428053357931682;
   0.59889391821413129 0.26240416483190626 0.5261024229156458;
   0.590115707522188 0.2564792889883718 0.51782061380025113;
   0.581197408128751 0.25062755963504108 0.50943804649957281;
   0.57213668952412411 0.24483471976895688 0.50096910722807086;
   0.562936972949324 0.23913740500288688 0.4924033487688122;
   0.55358965768684876 0.23348258786811352 0.48375741838397773;
   0.544127236829846 0.22794851634851182 0.47504601236622146;
   0.53453671498578681 0.22245473767373525 0.46623392778752137;
   0.52482925319759088 0.21705633258936946 0.4573608303884687;
   0.51500987065256942 0.21174067802933996 0.44842562979689171;
   0.50508369111093565 0.20651042518731044 0.439419706155005;
   0.49506769889479318 0.20131205349550388 0.43036462180177881;
   0.4849466656937394 0.19627886366835723 0.42124616339072041;
   0.4747642682051455 0.19127514845252946 0.41210007740167148;
   0.464496350234098 0.18639469150795496 0.40290052895775758;
   0.454153386568872 0.18157235929309856 0.39366805528979149;
   0.4437598142194083 0.17687622233814093 0.38441025973164494;
   0.43330800039241141 0.17225320322385185 0.37513214251476074;
   0.42282095722680396 0.16773266926121835 0.3658478935027954;
   0.41232198778425255 0.16332155241477755 0.35655335932520371;
   0.40177534359995631 0.15897329955709752 0.34725899080924172;
   0.39124856937586944 0.15471480093019349 0.33795541797906486;
   0.38071043058805748 0.15057593034381198 0.32869446559992754;
   0.37017496792619348 0.14651155490474935 0.31944634595997107;
   0.359685417924866 0.14257852861884387 0.31024521699190788;
   0.34923062289081058 0.13872067358942394 0.30105981157546496;
   0.33882932025805124 0.13498941162883121 0.29196431805269973;
   0.32849415810731691 0.13133013682463693 0.28292969271475782;
   0.31824140509122473 0.12777680038434783 0.2739626954515883;
   0.30808360186306172 0.12431088599079429 0.2650781883776554;
   0.29804703206605904 0.12096907737943005 0.25631109635573046;
   0.28815063647978806 0.11777548837573119 0.24768273091555792;
   0.2784138854051475 0.11462467691943093 0.23915689562797576;
   0.26884545195104881 0.11168663140944599 0.23078784280029102;
   0.25945683960653165 0.10877469369618373 0.222587654511519;
   0.2502458594990678 0.10605138371387143 0.21454770601388312;
   0.24130575993797584 0.10341309888075489 0.20672917281971837;
   0.23258145359967575 0.10085702750485699 0.19904990086445895;
   0.22409951967792754 0.098494129518777 0.1916317223679;
   0.21592800997000683 0.0961824842233464 0.18443431459626436;
   0.20798719274014368 0.094097878277804581 0.17748102188248838;
   0.20032199197668771 0.092102281739737252 0.17072404251303147;
   0.1929930519112138 0.090210148501606108 0.16425478778196315;
   0.18595756756550741 0.088461144321407864 0.15798916442765526;
   0.17918239607310721 0.086860604956501233 0.15197171824095515;
   0.17271652286342171 0.085310019059756326 0.14623268208498683;
   0.16657674927897742 0.084017159643581932 0.14075158682475447;
   0.16070109839239655 0.082744739570959444 0.13546266695617087;
   0.15514547215380362 0.081682831678010923 0.13049344769769583;
   0.14990229305804903 0.080652947861220028 0.12570362478180497;
   0.14492643188929413 0.079780198005186029 0.12112093754443654;
   0.14020433887741138 0.0790373210847673 0.11684795027326232;
   0.13577808508877007 0.078426472625953911 0.1128150277494095;
   0.13168472051183228 0.077943602177966 0.10894005513142714;
   0.12782469451000153 0.077585564996749457 0.10529427360161833;
   0.12421629669745488 0.07733153487505906 0.10190263239967523;
   0.12091263009869369 0.077160806294253037 0.09872438941971598;
   0.11792504160967056 0.077087605471277382 0.095739340788036892;
   0.11512428902632266 0.077124106477310786 0.092921438390002553;
   0.11266986288148367 0.077278280931717222 0.09034425711463924;
   0.11042054133233427 0.077557461916173034 0.087858259817374884;
   0.10835482664895428 0.0779682480970983 0.085430749270030623;
   0.10664980452486861 0.078515638239463337 0.083233362743108941;
   0.10500076149891276 0.079206706936435931 0.0811849987752464;
   0.10367568576805963 0.080047992441620724 0.079202373132585865;
   0.10245440459721128 0.081036146367393036 0.077408105510809136;
   0.10143302764234208 0.082172651362223986 0.075793213181368876;
   0.10060156744370947 0.0833429989095688 0.074343736373143038;
   0.099957490073521879 0.084733219859538786 0.073021493928374509;
   0.099491765663190221 0.086174200957104735 0.071799048152462019;
   0.099203561101382551 0.087868351782728449 0.070715865301635428;
   0.0990917059637016 0.089631418554732931 0.069812590208987113;
   0.0991543652151038 0.091581848333658483 0.069047413739548341;
   0.09938449662079421 0.093596705856330023 0.0683366602518671;
   0.0997592539364797 0.095870677385009817 0.067775575826204651;
   0.10029090453132976 0.098367735532022327 0.067351420714498469;
   0.10098624045100855 0.10100470443742694 0.067056376586652811;
   0.10185037074254032 0.10390316684470453 0.066890683170048468;
   0.10289669636944589 0.10701529424307119 0.066853405228907661;
   0.10407107749830058 0.11031323030702488 0.066942430823027158;
   0.10542706439764514 0.11379915823974537 0.067154572866541862;
   0.10700779523195278 0.11750256613139791 0.067485043150801907;
   0.10866356171507491 0.1214186354077631 0.06792871721692266;
   0.11058985460973997 0.12561122475299133 0.068490256434907332;
   0.11265358046302296 0.12998201016809363 0.069161658162312137;
   0.11482592212988228 0.13452952913401728 0.069841869202726614;
   0.1172525468096339 0.13923362540088058 0.07061013337583312;
   0.11984662822996572 0.14421992458906541 0.071527788775662821;
   0.12259023283587472 0.14936638834586771 0.07240250671869497;
   0.1255783308659662 0.15467097175200983 0.073463025646658617;
   0.12866632838600495 0.16015301428107631 0.074428953196786857;
   0.131964787295349 0.16583809622025805 0.075451143903149609;
   0.13539659003490429 0.17168797457807797 0.076499352947652935;
   0.13898143829032089 0.17771373035722049 0.077615266204072511;
   0.14273333408977723 0.18382162202537478 0.078814131582933572;
   0.14657648973463874 0.19010334825969472 0.0800978422742981;
   0.15058065647312086 0.19653643591917344 0.081472858980824625;
   0.1546789481267383 0.20304020599527237 0.082820325477988468;
   0.15891435033505735 0.20967894217083144 0.084315288412827941;
   0.16324201574954819 0.21644027415724562 0.085725671083784039;
   0.16764331569528698 0.22326088564616892 0.087378111733006275;
   0.17214117878407809 0.23014588067620151 0.088954666170349514;
   0.17672854171643612 0.23717030968648561 0.090616831164978665;
   0.18138557725174315 0.24418377345901152 0.092313787591837712;
   0.18614932716446331 0.25131965837076831 0.094070522565224732;
   0.19092496286475064 0.25845900126236249 0.095838932663733417;
   0.19578356565065036 0.26566646633545055 0.097701740737874521;
   0.20066624095575225 0.27289615431245268 0.099538585154345632;
   0.20563865734175896 0.28015543401981796 0.10144137335491779;
   0.21062165767952026 0.28743949794449813 0.10341681154350842;
   0.21565002131034081 0.29474849093360972 0.105341063504266;
   0.22071823626978779 0.30206589438149267 0.10737177098356988;
   0.22578580621308286 0.30942160308437588 0.10942425918828295;
   0.23087358793320908 0.31674637375346432 0.11146080916737203;
   0.23599942554298955 0.32407463157645766 0.113544148528805;
   0.24111667829948791 0.33140451979786362 0.11562706145018159;
   0.24625136795429642 0.33874382919190027 0.11774399033448354;
   0.25142161130991514 0.34605225161558245 0.11987530519597368;
   0.25655614356251627 0.35336798331579683 0.12201990271987771;
   0.26171486063228439 0.3606504869520773 0.12421841688792686;
   0.26686082811499112 0.36793415583869693 0.12644549680491185;
   0.27199976287771643 0.37518866036632931 0.12865270069350918;
   0.27716564862495718 0.38242312592009547 0.13091773299641;
   0.28230852522219579 0.38963517205482306 0.13315567648504587;
   0.28741102443064354 0.39681617951607306 0.13541131263839576;
   0.29253083679368164 0.40398253661323291 0.13772772056047208;
   0.29763371271500733 0.41110968925391994 0.13997703951686355;
   0.30270861597441245 0.41819983530612287 0.1423191896717482;
   0.30777796165744176 0.42527385370829907 0.14466141912573993;
   0.31283390340664441 0.4323092721765624 0.14698675045076079;
   0.3178717861969651 0.4392907679555289 0.14937445440835861;
   0.32289354884426558 0.44625078135606427 0.15173348027714889;
   0.32787458055379748 0.45317856550927788 0.15414019857003269;
   0.33286482215875229 0.46006214141313262 0.15659812478090693;
   0.3378065372786504 0.4669272611058824 0.15904172656258012;
   0.34276011776289156 0.47374251358940345 0.1615483986901225;
   0.34768618723430689 0.48054449680792832 0.16406525735832531;
   0.35260252743648796 0.48733177067501221 0.16660678833372614;
   0.3575330157602013 0.494096687049402 0.16922538605562049;
   0.36244695603233279 0.50086092275908345 0.17184748467098276;
   0.36737973646822791 0.50763595068013079 0.17458445512115095;
   0.37234266635346547 0.514431817878248 0.17737672866131532;
   0.3773460407454135 0.52125117454072178 0.18021855080909088;
   0.38237589730521832 0.52812192326384988 0.18317582716062439;
   0.38746281261163634 0.53504688383682231 0.18625771458158069;
   0.39261167353963283 0.54203578101426164 0.18942331524372824;
   0.397829085784384 0.5491052082225073 0.19271757689962568;
   0.40310886490283537 0.5562369035502045 0.19615957003504497;
   0.40846245878668191 0.56347569448539392 0.19969937055533532;
   0.41389916257600562 0.57078159533133987 0.20345060252485808;
   0.41942196889792105 0.57818748243299578 0.2073411155847612;
   0.42502574977129265 0.58569887098421314 0.211404841779217;
   0.4307139752414279 0.59329211638962209 0.21564607287093532;
   0.43649416325565288 0.60098043534195245 0.22009256839896127;
   0.44236981168434608 0.60878025364983912 0.22469947837233453;
   0.4483260263550442 0.61666512424328523 0.22956446645586148;
   0.45438943410546523 0.62465042027402884 0.2346800005644202;
   0.46052940028911515 0.632737035704805 0.23997201192980577;
   0.46679056980794126 0.640919568469155 0.24552927016246076;
   0.47312687628461886 0.64921252303499533 0.251375372684066;
   0.47958575169404788 0.65759918549801 0.257447237704721;
   0.48612398077214747 0.66607755923715506 0.263824123418753;
   0.49276575094068081 0.67466244733811409 0.27047165759239888;
   0.49951169685644409 0.68334645876598477 0.27740479242451704;
   0.50635983780977734 0.69213246608311907 0.2846367524995696;
   0.51331140915574447 0.7010071447066677 0.2922027961846701;
   0.52034930423735792 0.70998318036288943 0.30007680415059135;
   0.52749682832248757 0.719052417465989 0.30827939822186506;
   0.53473596274962532 0.728214506452531 0.31681629031981406;
   0.54206684387307391 0.73747096484237229 0.32567224953147406;
   0.54949627052740291 0.74681853385336983 0.33490548657583885;
   0.5570155260660562 0.75625216761439062 0.34443479006749489;
   0.56460896498657653 0.76577422583677857 0.35433684675399396;
   0.57229852335789722 0.775371785423672 0.36456981518560166;
   0.5800637024251929 0.78505863073597137 0.37515087340402431;
   0.58789315677452147 0.79481631444769063 0.38607083636919259;
   0.59580599779183585 0.80465115757879968 0.39733793110431614;
   0.60378647919259265 0.81455295139930606 0.40893684160743943;
   0.61181547035351136 0.82452569907804463 0.42085546058828815;
   0.61991079863233145 0.83456725158916312 0.43310537365549256;
   0.62804801875509753 0.84466593686424818 0.44565690798356666;
   0.63623203590160815 0.85482424353285191 0.45852386626204589;
   0.64444607574281032 0.86503411154502 0.47168426742096564;
   0.65270122774871442 0.875304603095295 0.48510702563508995;
   0.66098657874524047 0.88562373850225762 0.49881957040426217;
   0.66929937170304876 0.89599202339631379 0.5127757387973666;
   0.67763024722082377 0.90640898430139261 0.526986264888457;
   0.68597038266850219 0.91687395034085084 0.5414058019477237;
   0.694320711734304 0.927379717351087 0.55605017988029948;
   0.70269465411897525 0.93794159023020907 0.57089837356301387;
   0.71107079380329474 0.94854516351462181 0.58593239680412879;
   0.71944875140374809 0.9591979101370891 0.60111633537832909;
   0.72782482593864328 0.96989267598952911 0.61645627860765984;
   0.73620017000850424 0.98063388302955679 0.6319142852208457;
   0.74457595641659269 0.99141344595707814 0.6474770861182565];

% Default granularity
if nargin == 0
    N = 256;
end

% Rescale if necessary
if N ~= 256
    cmap = imresize(cmap, [N 3]);
end

% Ensure within range
cmap(cmap > 1) = 1;
cmap(cmap < 0) = 0;