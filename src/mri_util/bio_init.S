

BioSetFdes _ function(InFdes,OutFdes) {
  if( length(InFdes)>1 || length(OutFdes)>1 )
  stop("InFdes and OutFdes must be single integers")
  .C("bio_set_pipe_fdes",as.integer(InFdes),as.integer(OutFdes))
  1
}

BioReadFloats _ function(InBuf) {
  z _ .C("bio_read_floats",InBuf=as.single(InBuf),as.integer(length(InBuf)))
  z$InBuf
}

BioWriteFloats _ function(OutBuf) {
  z _ .C("bio_write_floats",OutBuf=as.single(OutBuf),
	as.integer(length(OutBuf)))
  z$OutBuf
}

