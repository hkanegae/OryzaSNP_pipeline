## OryzaSNP Pipeline をMacで実行する環境構築
Set up OryzaSNP Pipeline.

### Bioconda 設定
http://bioconda.github.io/
- conda がインストールされていない場合は、minicondaを入れる
- https://docs.conda.io/en/latest/miniconda.html
- 今回はMiniconda3を入れるので、Macの場合はMiniconda3 MacOSX 64-bit bashをダウンロードする
- ダウンロードしたシェルスクリプトを、ターミナルから実行
- sh ~/Downloads/Miniconda3-latest-MacOSX-x86_64.sh
- 設定は以下のページを参照してください
- http://imamachi-n.hatenablog.com/entry/2017/01/14/212719
- http://bonohu.jp/blog/bioconda.html



***

### vcftools のインストール
- biocondaを利用してインストール
- https://bioconda.github.io/recipes/vcftools/README.html
- conda install vcftools
- conda update vcftools

### vcftools　の実行
- parameter setting
  - max-missing=0.8 
  - max_alleles=2
  - min_alleles=2
  - maf=0.025
  - maxDP=1000
  - minDP=3


- vcftools --gzvcf sample.vcf.gz  --maf "$maf" --max-alleles "$max_alleles" --maxDP "$maxDP" --min-alleles "$min_alleles" --minDP "$minDP" --minQ 20 --max-missing "$max-missing" --out sample.fil --recode
- java -Xmx256g -jar 'beagle.28Sep18.793.jar'  gt=sample.fil.recode.vcf out=sample.fil.imputed



