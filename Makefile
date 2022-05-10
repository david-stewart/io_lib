all:
	./tu_make.sh

clear:
	rm src/*.pcm || true
	rm src/*.d   || true
	rm src/*ACLiC*   || true
	rm obj/*     || true
	rm src/*.so  || true
	rm lib/*     || true
	rm include/* || true

