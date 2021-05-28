all:
	./make.sh

clear:
	rm src/*.pcm || true
	rm obj/*        || true
	rm src/*.so  || true
	rm src/*.d   || true
	rm lib/*        || true
	rm include/*    || true
