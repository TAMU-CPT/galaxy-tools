curl \
	-H 'Connection: close' \
	-H "REMOTE_USER: $1" \
	http://10.42.170.10:8999/apollo/user/checkLogin
