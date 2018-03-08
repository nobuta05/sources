# 作業用ディレクトリ
wd=$(mktemp -d);
# 作業用ファイル
temp=$(mktemp)

number=1
while read text; do
  if [ -n ${text} ]; then
    id=`printf "%04d" ${number}`
    echo ${text} > $wd/inp.txt;
    open_jtalk -m /usr/share/open-jtalk/voices/nitech_jp_atr503_m001.htsvoice -x /usr/share/open-jtalk/dic -ow $wd/${id}.wav $wd/inp.txt;
    lame -V2 $wd/${id}.wav $wd/${id}.mp3
    rm $wd/inp.txt $wd/${id}.wav;
    number=$(( ${number} + 1 ))
  fi
done < $1

sox $wd/*.mp3 $2"/"$3"_"${1%.*}".mp3"

##### 終了処理
# rm -rf $wd
rm $temp