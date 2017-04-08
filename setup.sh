wget http://meme-suite.org/meme-software/4.11.3/meme_4.11.3_1.tar.gz
tar -xvzf meme_4.11.3_1.tar.gz
cd meme_4.11.3
./configure --prefix=$HOME/meme --with-url=http://meme-suite.org --enable-build-libxml2 --enable-build-libxslt
make
make test
make install
export PATH=$HOME/meme/bin:$PATH 

