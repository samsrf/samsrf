function cmap = romaO(N)
% Colour map by Fabio Crameri: 
%   Crameri, F. (2018), Scientific colourmaps, Zenodo, doi:10.5281/zenodo.1243862
%

cmap = ...
  [0.45137387466634432 0.22345870991841682 0.34187079996534681;
   0.45417936034863055 0.22244423818189635 0.3360991540854924;
   0.45696452229328971 0.22158360152233506 0.33042757476772827;
   0.45974582046150508 0.22090162872323899 0.32483395049002461;
   0.46250875806408059 0.22034558424660483 0.31934641011548964;
   0.46526890288254408 0.21993574221232304 0.31393778773564368;
   0.46802588595732719 0.21968046217604173 0.30861938574464759;
   0.470782101675561 0.21957701023228535 0.30336671980434493;
   0.47351926572639552 0.21962388229695 0.29821777373013758;
   0.47627535151613221 0.21982136967640389 0.29315815972377174;
   0.47902403765590024 0.22017173928147912 0.28818364522263229;
   0.48177891799147871 0.22067262231068385 0.28330406637410249;
   0.484529958777781 0.22130249026628343 0.27850233094524635;
   0.4873079402723392 0.22207902216818232 0.27379001398985814;
   0.49007871542098724 0.22304082989188881 0.26916688247332315;
   0.49286240296059319 0.2241064037394376 0.26460608433016053;
   0.49567010128681244 0.22536308195989443 0.26015913280103442;
   0.49850015244630486 0.22677002702174462 0.25578605964674661;
   0.501335001941077 0.22832529695482012 0.2515283329633044;
   0.50419462152951655 0.22999167531668241 0.2473312990543006;
   0.50706548581322941 0.23188016281160506 0.24321985175414734;
   0.50996795728322664 0.23386578649577738 0.23922939196602713;
   0.51289632809273045 0.23604989933815851 0.23533362221163814;
   0.51584420069121761 0.23835010381262256 0.23151439528315168;
   0.51884191184459216 0.24082274406310569 0.22778705171241884;
   0.52184159401472685 0.24345287921714415 0.22413643277264281;
   0.5248936633081952 0.24624614345720136 0.22064654139407122;
   0.52796880060702112 0.24920033871202815 0.21720293269271293;
   0.53108222657902271 0.25230697241545824 0.2138733240350868;
   0.53422500404267692 0.2555629709908836 0.21064330259694405;
   0.537416087128431 0.25898718755867778 0.20752821570454091;
   0.54062953663424618 0.26254874906527853 0.20452422558065472;
   0.54388578501156071 0.26628074194001894 0.20158479795845374;
   0.547179477632665 0.27016871338346554 0.19879103455372155;
   0.55051470417601389 0.27419220768413571 0.19613158941587366;
   0.55388507057000047 0.27838557801396008 0.19355896640738146;
   0.55730513843560336 0.28273496624487637 0.19108901259165578;
   0.5607481891849585 0.28720295365543952 0.18877025772830527;
   0.56424420538805009 0.29186044607623146 0.18654823602603374;
   0.56777225640943729 0.29664943171883 0.1844613681005102;
   0.5713395708372605 0.30156762515375779 0.18248255558096546;
   0.57494760412139745 0.306664119697607 0.1806547327036454;
   0.57859849037404365 0.3118638725943903 0.17897706912758035;
   0.58228385364368362 0.31724005709857195 0.17742974625817579;
   0.58601668474934554 0.32274917620785049 0.17596600443205185;
   0.58977086382260246 0.32837579170005987 0.17473041913499696;
   0.59357770099990326 0.3341520737807156 0.17358461186773633;
   0.59741621424355873 0.34005082962595734 0.17260771790913279;
   0.60128520864126933 0.34606488732098356 0.17179107102537722;
   0.60519429347073328 0.35222612613486654 0.17113882695994312;
   0.60914595966965834 0.35850660092723158 0.1706518001661704;
   0.61311403585734081 0.36491230545545161 0.17033659138461416;
   0.61713305263738794 0.37143170001699405 0.17019647243974415;
   0.62117831817168534 0.37808004260512812 0.17023424723468036;
   0.62525650405948174 0.38483277844679603 0.17045584576837433;
   0.62936663645221713 0.39171274010832574 0.17086949927169515;
   0.63351583534698053 0.39869161366526606 0.17148029832673492;
   0.6376899801592133 0.40579273451828446 0.17229455415902861;
   0.6418969206640357 0.41299341858588462 0.1733217031052732;
   0.64613478038363159 0.42028968631183211 0.17457947831402124;
   0.65040731118022488 0.42770821256165048 0.1759997030861708;
   0.6547034678186211 0.43522198057514994 0.17773555701052837;
   0.65903594688945388 0.44283104701509224 0.17961872462744535;
   0.66340663709919723 0.45053919844812856 0.18175410244602871;
   0.66779754864849228 0.45834400246563994 0.18415883675091;
   0.67221701556636193 0.46624601288315842 0.18679653303009922;
   0.676668408258122 0.47424556516199806 0.18968343894008491;
   0.68113987453521319 0.48233432047961189 0.19282559174242878;
   0.68565504614384887 0.49051384176748442 0.19624212364825352;
   0.6901865943841935 0.49877525247161042 0.19987151086794824;
   0.69474489267519313 0.50711878939777111 0.20384095613897149;
   0.69933150190502635 0.51554364510315542 0.20802951717189666;
   0.70394379944458541 0.52405913308757035 0.21250815807665588;
   0.708579269984744 0.53264810851365207 0.21726466652424309;
   0.71322462487202609 0.54129584783370643 0.22228801954461327;
   0.71789543172411874 0.55002850287196436 0.22760999789642963;
   0.72257074879630612 0.55881191999363145 0.23318038744449734;
   0.72726541908143694 0.567667619090913 0.23906844733572025;
   0.73196754124822017 0.5765770306810214 0.24521274313140504;
   0.73666429855520665 0.58552606964382192 0.25167690199535026;
   0.74135966280553234 0.59451143602026024 0.25837336184024506;
   0.74605208722083227 0.60354190354601767 0.26536942297663035;
   0.75073206970193451 0.612586769528166 0.2726302004336314;
   0.75538458761417815 0.62165641951703632 0.28016536652975538;
   0.76000626375777725 0.63074519072858415 0.28796159285400263;
   0.76459850548382058 0.63982107372195873 0.29602133603704045;
   0.7691393547924017 0.6488872007879265 0.30432641333726018;
   0.77363082353747292 0.6579334501062053 0.31287317590766928;
   0.778063513547249 0.66693648173358677 0.32165342045778461;
   0.78241858181306168 0.67590452134349488 0.33065628764834415;
   0.78669473222742692 0.68481012848910849 0.33988429768103184;
   0.79086729978447268 0.69364601640968437 0.34929177005468748;
   0.79494152148162034 0.7024000734374567 0.35888066437276217;
   0.79890056583658753 0.71106110014509849 0.36866666904924328;
   0.80272657351953081 0.7196075866437891 0.37858785537462053;
   0.80642183739287643 0.7280278386545781 0.3886641443829395;
   0.809957948959543 0.73631442291646487 0.39885356439890424;
   0.81333508775589614 0.7444579033381834 0.40916419343044097;
   0.8165453805484858 0.75243576471807627 0.4195678589512597;
   0.81956066662524385 0.76024551955216324 0.43003909883623515;
   0.82239352145030631 0.76786527606303778 0.44056580600378381;
   0.82501007998665521 0.77529503145242462 0.45114566396597466;
   0.82741842101803253 0.78252062840137226 0.46174207144418739;
   0.829603946670109 0.78953125251921252 0.47234557018298379;
   0.83155421929555018 0.79631345904307094 0.48293113505348895;
   0.83325763948187437 0.80286945546535671 0.49349438092320991;
   0.83471711440316687 0.80919217337921789 0.50402349238016431;
   0.83591802992904807 0.81526290648804256 0.51448637001357589;
   0.83685834726673813 0.82108747825788586 0.52487063926499433;
   0.83753112877877889 0.82666152911953383 0.535169715855738;
   0.83792711883984428 0.83197775136575325 0.5453749940327054;
   0.8380452219661233 0.837030911276995 0.55545503373542693;
   0.837882936293263 0.84182372978083053 0.56541687164200261;
   0.83743505691465325 0.84634841066352917 0.57525009248389325;
   0.83669608738001033 0.85061428946982776 0.58493481683268178;
   0.83566900297218039 0.85461822150317668 0.59445662475083871;
   0.834353490535377 0.858352935600571 0.603824386219959;
   0.83273908958133624 0.86183106398244691 0.613005156127519;
   0.83083916784666578 0.86503881681958983 0.62201642584302752;
   0.82864273153821466 0.867996924203077 0.630853456571219;
   0.82615289276805282 0.87068453138380941 0.63948622544611522;
   0.82337316226109114 0.8731243808170851 0.64792166272266227;
   0.82030302049787251 0.875308817916059 0.65617051431273388;
   0.8169502875118102 0.87723785240565988 0.6642016059109499;
   0.81330617911779668 0.87892193607412561 0.67203050372942275;
   0.809388073505355 0.88035883911100588 0.67964369368348843;
   0.80518103676978658 0.88155539547271411 0.68705262496163488;
   0.80070557418553945 0.88250384232711687 0.69423702355165373;
   0.79595416986815615 0.88322024160776869 0.70121490825136612;
   0.79093896722814139 0.88370166442568776 0.70796955376492232;
   0.78566354355748846 0.88394756791192053 0.71450239681921568;
   0.78012448316467631 0.88396267845785781 0.720818601436446;
   0.7743349809549066 0.88374961537775332 0.72691514848868211;
   0.76830117500491169 0.88330880657064237 0.73279235943307086;
   0.76202932867127926 0.8826439664605501 0.73844080219926611;
   0.75552745305754587 0.88176528543097576 0.74387308004301733;
   0.7487922201442716 0.88066142855267426 0.74907901244010133;
   0.741840965356767 0.87934274252982469 0.75407362506942921;
   0.73468289435733047 0.877813133624419 0.758841121047384;
   0.72731286344152735 0.8760681327994988 0.7633884719064522;
   0.71975563131740072 0.87411281342913594 0.76771863702656429;
   0.71200634254652917 0.87195490205629766 0.77183731326807825;
   0.70408157778278535 0.86958361586747523 0.77573189149537913;
   0.6959934435841667 0.86701199542454632 0.77941408611244967;
   0.68774086401590184 0.86424665083092067 0.78288161414727053;
   0.6793389162237915 0.8612733805955386 0.786142515905421;
   0.67081043488060488 0.858108754869297 0.78918696963758916;
   0.66214974823204054 0.854755325633921 0.792019374606443;
   0.65336191615535322 0.8512038147137011 0.79464949462704715;
   0.64447888526236663 0.84747281183107692 0.79706929443474994;
   0.63550466648887127 0.84355688873942836 0.79929550137047922;
   0.62644789368993981 0.83946736643345121 0.80131005301653169;
   0.61732250245054976 0.83519465672573878 0.80313031334474994;
   0.60814458330216625 0.83075393097610029 0.80475747344017912;
   0.59891207904994292 0.8261419477718408 0.80618974025285306;
   0.58964582362531281 0.821370001503925 0.80743073201209858;
   0.58037416891668359 0.81644069929738694 0.80848448431178188;
   0.57108295795564334 0.81135417052093117 0.80935499467904015;
   0.56181096656633911 0.8061197142563109 0.81003830104829577;
   0.55254696050788421 0.80074153357692224 0.81054737180923975;
   0.54331692938890741 0.79522487473470926 0.81088486803854054;
   0.53412481343869844 0.78957864780153264 0.81105073612913037;
   0.52499942801378252 0.78380494339557993 0.81104829255541788;
   0.51592946640408521 0.77790619684681317 0.81088128535264992;
   0.50695457178679948 0.77189237606949457 0.810553867312021;
   0.49807730814939721 0.76577065591917992 0.81007313367329126;
   0.4892761154561634 0.75953627145039693 0.80944413096194612;
   0.48060724428672041 0.75321035455809271 0.80866209982050219;
   0.47207158472482136 0.74679509936917687 0.80773478782978692;
   0.4636543832157532 0.74029003314965958 0.80667103569997511;
   0.455394360620785 0.733700555344614 0.80546175528685726;
   0.44727892783103695 0.72703421775700483 0.80413015801056875;
   0.439336620501912 0.72029850538486406 0.80266092075400075;
   0.43157820295141064 0.71349722048094144 0.80107034066150862;
   0.42397607905257317 0.70663759445882113 0.79936048983979813;
   0.41658319040800768 0.699710836281273 0.79752058436818341;
   0.40937950998789241 0.692746942531493 0.79557028868455426;
   0.40237222654930827 0.6857226343304319 0.79351148257178772;
   0.39559768851882859 0.67865464677983212 0.79132980443516354;
   0.38902594754719722 0.67154890742151208 0.78904829753662376;
   0.38267134510866385 0.664406303318654 0.78665686031850712;
   0.37656194726874181 0.6572364599940409 0.784158256900166;
   0.37065998671735456 0.650026169816983 0.78155312432186663;
   0.3650156478374445 0.64278907262214235 0.7788439704346074;
   0.35961474357462914 0.63552459531347538 0.77603609574121668;
   0.35446128756534578 0.62824435983995219 0.77312335404412325;
   0.34955329783298811 0.62093578001801386 0.77011259356745276;
   0.34490132688387537 0.61359860829748269 0.76699764699443562;
   0.34051206148341034 0.60625349053160149 0.76378192120145139;
   0.33637407306721356 0.59889333087048446 0.76046535323702757;
   0.33252896897999634 0.59151385914039534 0.75704026602868879;
   0.32892978105966453 0.58411541639695075 0.753509660145864;
   0.32558660956763974 0.57671386227694665 0.7498713881606146;
   0.32255538145454354 0.56928040025995874 0.74612904377833811;
   0.319776812327851 0.5618552336234679 0.74227540162875638;
   0.31727215769941852 0.5544059687934092 0.73830238489976308;
   0.315046235440496 0.54694880751294206 0.73422056502553468;
   0.31310888100727757 0.53947683250337009 0.730018305046634;
   0.31144078510148449 0.532009191399311 0.72569495131483008;
   0.31007085573219206 0.524525177713516 0.72123940516733365;
   0.30896596754459188 0.51704235142424815 0.71667268924133032;
   0.30811233694038503 0.50954589542565287 0.71196550186717888;
   0.30754755596155547 0.5020535600125019 0.70713165248419552;
   0.30725992941660535 0.49456098903836859 0.7021584317490277;
   0.30723357676160024 0.48706527979628977 0.69705579801043982;
   0.30746278637107111 0.47957805386820157 0.69181631044117931;
   0.30794945480608615 0.47209780815860447 0.68643200203671917;
   0.30870214160456039 0.46462786076752843 0.68090429224240789;
   0.30967656288385531 0.45716140421572532 0.67524510365912682;
   0.3108794327886486 0.44972766840391859 0.66944262233617879;
   0.31227930455350716 0.44231871442218229 0.66349906912858048;
   0.31393417849302185 0.43492966501369112 0.65741193303837242;
   0.31577838324933488 0.42757684343728208 0.65118353596296263;
   0.31779510113080289 0.42025287664636152 0.64481659797988367;
   0.32000908101745906 0.4129924331824989 0.63832609199652368;
   0.32237721743409953 0.40577354024135681 0.63170368958752321;
   0.32488595189411384 0.39860600863872309 0.62495431234654519;
   0.3275540524524464 0.391515785986681 0.61808640021287031;
   0.33035290864531491 0.38448361014541277 0.61110958069392574;
   0.33326998516469519 0.37754527439182789 0.60402426329884329;
   0.336267009070308 0.370682136446979 0.596835463446233;
   0.33938746876761527 0.3639232716213518 0.58954910977947228;
   0.34256617055668848 0.35727961323829893 0.5821911715106195;
   0.34579964982296035 0.35072840114129222 0.57475556754020807;
   0.34911680950670931 0.34427606684739553 0.56726723235531129;
   0.35246712188286733 0.33796951086642585 0.55971421876734218;
   0.35586571330385003 0.33179429361522866 0.55211951999674158;
   0.3592719698405375 0.32574339895086335 0.544484280602014;
   0.362710136498347 0.319862449027072 0.5368357564177213;
   0.36616701523356 0.3141059686703076 0.52916547338239317;
   0.3696106729341827 0.30851510115161712 0.52148082645634475;
   0.37306065699548407 0.30306340153795336 0.513820682513332;
   0.37651721415497957 0.29779593012417804 0.50615476903713374;
   0.37993913271044283 0.29269060993945273 0.49853693659327875;
   0.3833601695112368 0.28775211045886473 0.49094098365969924;
   0.38674034893865988 0.28300857052541106 0.483366122680233;
   0.39011148602639573 0.27841787444716054 0.47586203774627034;
   0.39345972212767938 0.27400888505045051 0.46839655688640258;
   0.39677187524832774 0.2697844478818473 0.46099821935897295;
   0.40006389871098774 0.26572778671695885 0.4536648650121633;
   0.40332570708290255 0.26185063531041125 0.44640357258705887;
   0.40654740140736317 0.25814609463915134 0.43921282229253944;
   0.40974122584995359 0.25465794070620523 0.43211955213998066;
   0.41290484428531121 0.2513205392469372 0.42508691943295418;
   0.41602298069925253 0.24816571088231248 0.41813476334746141;
   0.41912163570789673 0.24515044941296285 0.41128263957791467;
   0.42217639304112342 0.24235173234917753 0.40451146776521518;
   0.42522374243151861 0.23972006660417478 0.39784449951144712;
   0.42822601985526054 0.2372759110590083 0.39126117987526732;
   0.43120703766852425 0.23497723471652363 0.38475458712696586;
   0.43414510340680063 0.23282014437738366 0.37835660182672837;
   0.43708174492247281 0.23086408171643055 0.37204042423124;
   0.4399789158880909 0.22906557783149017 0.36583161097363492;
   0.44285961094407655 0.22742791946705074 0.35970303621083011;
   0.44570569896935458 0.22595659601399223 0.35366131629413211;
   0.44855465158106778 0.22459614104497652 0.34772837453617844];

% Rotate by 90 deg
cmap = [cmap(65:256,:); cmap(1:64,:)];
% Flip vertically
cmap = flipud(cmap);

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