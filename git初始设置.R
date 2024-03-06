
#设置路径
cd "C:/Users/21989/Documents/test.A"
git init
git config --global --list
git remote -v
git remote add origin https://github.com/521LWS/test.A.git
git push -u origin master # 同git push --set-upstream origin main，相当于设置默认推送分支
#设置账户
git config --global user.email "2628220830@qq.com"
git config --global user.name "521LWS"
#设置密码
gpg --full-generate-key #生成密钥
gpg --list-keys #查看密钥
gpg --armor --export 33DBB18B01C92DD677C2AF74AD886BCCCE17DC51 #生成符合github格式的密钥
git config --global user.signingkey 33DBB18B01C92DD677C2AF74AD886BCCCE17DC51
git config --global commit.gpgsign true
git config  --global credential.helper manager
#定义pull后merge方式
git config --global merge.default recursive
git config --global push.default simple #多余但保险
.................................................................................
#处理冲突
git config --global push.default simple #多余但保险
git push --set-upstream origin main
git stash
git pull origin main --allow-unrelated-histories
git stash pop
git add git初始设置.R

...............................................................................
#删除历史记录
#切换分支：
git checkout --orphan latest_branch
#添加到暂存区：
git add -A
#提交更改： 
git commit -am "commit message"
#删除分支： 
git branch -D main
#重命名分支： 
git branch -m main
#强制提交到远程仓库：
git push -f origin main 
#重新连接跟踪分支
git branch -u origin/main
git stash
git pull origin main --allow-unrelated-histories
git stash pop
...........................................................................
#删除历史记录
#切换分支：
git checkout --orphan latest_branch
#添加到暂存区：
git add -A
#提交更改： 
git commit -am "commit message"
#删除分支： 
git branch -D main
#重命名分支： 
git branch -m main
#强制提交到远程仓库：
git push -f origin main 
#重新连接跟踪分支
git branch -u origin/main
#重新基于远程分支
git rebase origin/main
#终止 rebase：
git rebase --abort
#不管历史分支重新pull
git pull origin main --allow-unrelated-histories
git push -u origin main
git branch -D master #删除本地分支 （D 强制）
.......................................................
#唤醒自动merge
git pull origin main --allow-unrelated-histories
git push -u origin main
git push -u origin main
git push -u origin main
git push -u origin main
git push -u origin main
git push -u origin main
git config --global --list


