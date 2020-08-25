echo "# python" >> README.md
git init
git add README.md
git commit -m "first commit"
git config --local user.email "1873350971@qq.com"
git config --local user.name "jiarf"
git config --list
#git config --unset user.email
#git commit -m "first commit"
#git config --local user.email "1873350971@qq.com"
#git config --local user.name "jiarf"
#git config --list
git remote add origin https://github.com/jiarf/python.git
git push -u origin master

git add temp.py
git commit -m "Add temp.py"
git push origin master
#输入名字和密码（常用）

##若有一个项目很大，都在一个文件夹里，数据太大上G了导致无法上传，此时需要建一个.gitignore
vim .gitignore
test3.py
q


#在命令行中进行
git commit -m "test"
#时test3.py就会被忽略，不会被添加




#ignore everything
*
！.gitignore
!*.py
!*.pl
!*.R
!*.README
!*.sh
!*.ipynb
#
!*/
tool/


###
$ git add .
$git commit -m "first back"
$git push origin master
$git remote -v
####找到fetch的链接，liu览器打开

##$history | grep 'echo'

#readme.md中可以显示出来#一级目录##二级目录更改后
$git add README.md
$git commit -m "updata README.md"
$git push origin master














