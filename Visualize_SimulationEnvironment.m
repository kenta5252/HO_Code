%% 理論解析元論文_Vertical Handover Analysis for Randomly Deployed Small Cells in Heterogeneous Networks
%% Simulation のプログラミング

%% 円の中に基地局をランダムに配置

clc;
clear;tic;

%% パラメータ値設定
K=10000; %最大発生乱数の個数
NR=625; %観測エリア半径R内に発生させたいSセルの数(乱数の個数)
In_S=0; %円の中の基地局数(乱数)
Out_S=0; %円の外の基地局数(乱数)
R=1.0; %観測エリアの半径[km]
r=0.04; %Sセル半径[km]

%% 円の中に基地局をランダムに配置するための準備
numbSim=2; %x,y座標で保存するため
lambda=NR/(pi*R^2); %ポアソン点プロセスの強度λ(平均密度)=発生させたい個数/面積
xMin=-R; xMax=R; yMin=-R; yMax=R; %観測エリアの寸法
xDelta=xMax-xMin; yDelta=yMax-yMin; %差
% 一度にすべての設定を生成する
numbPointsA=poissrnd(NR,numbSim,1); %点のポアソン数 numbSim行1列のポアソン発生
numbPointsTotal=sum(numbPointsA); %ポアソン値の合計 ex.1221
x=xDelta*(rand(numbPointsTotal,1))+xMin;%ポアソン点のx座標 xDelta、yDeltaで規格化 xMin、yMinでオフセット
y=yDelta*(rand(numbPointsTotal,1))+yMin;%ポアソン点のy座標 生成されたランダムな座標を観測エリア内に収め観測エリアの左下隅を原点とする座標系に変換するための処理

%円の中の乱数を探す
%2×2の中にポアソン過程に従ってランダム配置をさせることで、とにかく発生させて、NRの数に達したらもし円の中にSセルがあっても無視
%　→　別に円の中にNR個のSセルがポアソン過程に従っているので問題ない
%　でもこのコードだとNR×2個を2×2の中に発生させようとしちるので、本当にNR×2であるひつっ用があるのかは考える必要あり
for N_in=1:K
    if ((x(N_in)/R)^2+(y(N_in)/R)^2)^(0.5)<=1
        ransu_naka(In_S+1,1)=x(N_in);
        ransu_naka(In_S+1,2)=y(N_in);
        In_S=In_S+1;
    end 
    if In_S==NR
        break
    end
end
%円の外の乱数を探す
for N_out=1:N_in %N_in(円の中の乱数を数えたのと同じ数)まで
    if ((x(N_out)/R)^2+(y(N_out)/R)^2)^(0.5)>1
        ransu_soto(Out_S+1,1)=x(N_out);
        ransu_soto(Out_S+1,2)=y(N_out);
        Out_S=Out_S+1;
    end    
end
%ransu_naka
%ransu_soto
%In %観測エリア内(円内)
%Out %観測エリア外(円外)
%% グラフ作成
%マクロセルの円作成
angl2rad = pi/180;
dataN = 100;
theta=linspace(0,360,dataN).*angl2rad;  %y=linspace(x1,x2,n) は、x1~x2 の間の等間隔の点をn個含む行ベクトルを返す
h=figure(1); clf;
x_centre = 0;  y_centre = 0;
x=R.*cos(theta) + x_centre;
y=R.*sin(theta) + y_centre;
plot(x,y,'red');    
%スモール基地局(乱数)をプロット
hold on
for n=1:In_S
x=ransu_naka(n,1);
y=ransu_naka(n,2);
scatter(x,y,'blue','.')
end
for n=1:Out_S
x=ransu_soto(n,1);
y=ransu_soto(n,2);
scatter(x,y,'cyan','.')
end
%観測エリア内スモールセル作成(円作成)
angl2rad = pi/180;
dataN = 100;
theta=linspace(0,360,dataN).*angl2rad;
for n=1:In_S
x_centre=ransu_naka(n,1);  y_centre=ransu_naka(n,2);
x=r.*cos(theta) + x_centre;
y=r.*sin(theta) + y_centre;
plot(x,y,'blue');  
end
%観測エリア外スモールセル作成(円作成)
angl2rad = pi/180;
dataN = 100;
theta=linspace(0,360,dataN).*angl2rad;
for n=1:Out_S
x_centre=ransu_soto(n,1);  y_centre=ransu_soto(n,2);
x=r.*cos(theta) + x_centre;
y=r.*sin(theta) + y_centre;
plot(x,y,'cyan');  
end
%マクロセルプロット
x=0; y=0;
scatter(x,y,'red','*')
%hold off
title('想定環境')
xlabel('x[km]') 
ylabel('y[km]')
%% UEが移動する

%1単位で動く距離を決める(細かさ)
tani_ugoki=0.001; %0.0001で十分かな
%始点U→終点Vのペアをランダムに決める
kakudo_U=2*pi*rand(); %始点のランダムな角度
kakudo_V=2*pi*rand(); %終点のランダムの角度
%極座標でまず0～2πでランダムに決めた後に、直交座標に直す
x_U=R*cos(kakudo_U); y_U=R*sin(kakudo_U); %始点の直交座標
x_V=R*cos(kakudo_V); y_V=R*sin(kakudo_V); %終点の直交座標
scatter(x_U,y_U,'yellow','p');
scatter(x_V,y_V,'yellow','p');

S=((x_V-x_U)^2+(y_V-y_U)^2)^(1/2); %U-V間距離
Cos_sita=(x_V-x_U)/S;
Sin_sita=(y_V-y_U)/S;

for n=1:In_S
kaunto(n)=0; %=1なら過去に接続済み、=0なら過去に接続なし
end
Ho=0; %HO回数
Vho=0; %VHO回数
E_N=0; %sセル通過数
joutai=2; %今の状態=2ならmセル、3ならsセル、初期状態はmセル
x_now(1)= x_U ; %ループのために初期状態を配列にいれる
y_now(1)= y_U ; %ループのために初期状態を配列にいれる
for nn=2:1000000002 %nn回単位当たりで動く
x_now(nn)= x_now(nn-1) + tani_ugoki*Cos_sita;
y_now(nn)= y_now(nn-1) + tani_ugoki*Sin_sita;
if nn>(S/tani_ugoki) %このifでnnのforループを終わらせるように
    break;
end 
    for n=1:In_S %スモールセルの数In個をnで全てみる
        Length_scell_user=((ransu_naka(n,1)-x_now(nn))^2+(ransu_naka(n,2)-y_now(nn))^2)^(1/2);
                    if Length_scell_user <= r
                        if kaunto(n)==0 %そのセル内に過去滞在ないならばHOする
                        Ho=Ho+1;
                        kaunto(n)=1;
                        E_N=E_N+1; %sセル通過数のカウント
                            if joutai==2 %1つ過去時間がmセルならば
                                Vho=Vho+1;
                            end
                        joutai=3; %状態をsセルにする    
                        break;
                        else %過去滞在あるならばHOしない
                            break;
                        end
                    end       
    if n==In_S %今現在の位置にsセルがないとき(Inすべて見た)→mセルにつないでい行く
      if joutai==3 %過去がsセルならば
          Ho=Ho+1;
          Vho=Vho+1; %Vho回数を増やす
          joutai=2; %状態をmにする
          break;
      else %過去がmセルならば
        break;
      end
    end
    end
end

Ho  %今回は参考程度(正確なものは求めてない、sセル→sセルのHOをいつやるのか決めてない)
Vho
E_N
In_S
%% ユーザの通り道を可視化(プログラミングをforで何回も回す時は省く。)
%hold on
for nn=1:1000000002 %nn回単位当たりで動く
if nn>(S/tani_ugoki) %このifでnnのforループを終わらせるように
    break;
end
scatter(x_now(nn),y_now(nn),'green','o');
end
hold off
%以上のプログラミングをforで何回か回して平均を取る


toc;