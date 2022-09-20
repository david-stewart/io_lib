tu:
	./tu_make.sh
all:
	./make.sh
_class: 
	./tmake.sh

clear:
	rm src/*.pcm || true
	rm src/*.d   || true
	rm src/*ACLiC*   || true
	rm obj/*     || true
	rm src/*.so  || true
	rm lib/*     || true
	rm include/* || true

