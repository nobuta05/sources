# 1. pdfファイルならページ毎に分解. pdf -> jpg
# 2. jpgをOCR
# 3. ページ毎にtxtを音声データ化

# 作業用ディレクトリ
workdir=$(mktemp -d);
# 作業用ファイル
tmpfile=$(mktemp)

# ファイルかどうかを判別
# pdfファイルなら分割してjpgに
# ディレクトリなら画像ファイルを作業ディレクトリに移す.
if [ -f $1 ]; then
  if [ ${1##*.} = pdf ]; then
    filename=${${1##*/}%.*};
    pdfseparate $1 ${workdir}/%03d.pdf;
    for f in $(ls ${workdir}/*pdf); do
      convert ${f} ${f%.*}".jpg";
    done
  else
    exit 0;
  fi
elif [ -d $1 ]; then
  filename=${1##*/};
  cp ${1}/*.jpg ${workdir};
fi

# 画像ファイルの余白処理
for f in $(ls ${workdir}/*.jpg); do
  mogrify -fuzz 50% -trim ${f}
done
julia /home/nobuta05/Gits/sources/template.jl ${workdir};

# OCR処理
for f in $(ls ${workdir}/*.jpg); do
  id=`gdrive import ${f} | awk '{print $2}'`;
  gdrive export --mime text/plain ${id};
  mv ${f##*/}".txt" ${workdir}
  gdrive delete ${id};

  nkf --overwrite -Z ${f}".txt"
  sed -i '/[。|、|0-9]/!d' ${f}".txt"

  echo "Finish "${f}
done

# 章ごとに分割
alltxt=$(mktemp ${workdir}/XXXX.all);
# ls $workdir | grep ".txt" | xargs cat -s > $alltxt
cat -s ${workdir}"/"*".txt" > ${alltxt}; # 1つ上のコマンドは，このコマンドの代わり
sed -i -r 's/ //g' $alltxt;
sed -i -r 's/\[([0-9]+)\]/\n\[\1\]/g' $alltxt;
sed -i '/^$/d' $alltxt;
# rm $workdir/*txt

linens=($(grep -n -E '\[[0-9]+\]' $alltxt | awk 'BEGIN{FS=":"}{print $1}' | xargs))
num=1
for i in $(seq ${#linens}); do
  # 最終行の処理
  if [ $i = ${#linens} ]; then
    tail -n +${linens[$i]} $alltxt > $workdir/$(printf "%03d" ${num})".txt"
  # 最終行以外の処理
  else
    # 第一要素が1でない場合の処理
    if [ $i = 1 ] && [ ! ${linens[$i]} = 1 ]; then
      head -n $((${linens[$i]} - 1)) $alltxt > $workdir/$(printf "%03d" ${num})".txt"
    else
      sed -n ${linens[$i]},$((${linens[$(($i + 1))]} - 1))p $alltxt > $workdir/$(printf "%03d" ${num})".txt"
    fi
  fi
  num=$(($num + 1));
done

# 章毎のファイルを音声データに変換
for f in $(ls $workdir/*txt); do
  # 各ファイルを一旦一行にまとめて，読点ごとに改行
  cat $f | tr -d '\r\n' > $tmpfile
  cat $tmpfile > $f
  sed -i 's/。/。\n/g' $f

  # 各ファイルを音声データに変換
  source /home/nobuta05/Gits/sources/generate.sh $f $workdir $filename
done

##### 修了処理
# rm -rf $workdir;
rm $tmpfile