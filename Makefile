all:
	./make.sh

clear:
	rm src_io/*.pcm || true
	rm obj/*        || true
	rm src_io/*.so  || true
	rm src_io/*.d   || true
	rm lib/*        || true
	rm include/*        || true
