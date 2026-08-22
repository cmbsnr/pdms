for file in *.data
do
	name=${file%.data}
	echo $name
done
