---
title: "使用hugo建立网站"
date: 2020-04-05T08:58:51+08:00
draft: false
---

> hugo轻量化，无依赖，速度快的静态博客生成程序。
> 个人博客预览：[https://xinli0111.gitee.io](https://xinli0111.gitee.io)


# Hugo下载
Hugo官网：[https://gohugo.io/](https://gohugo.io/)<br />安装方式参考官网

# 开始建立网站

## 1. 建立网站
建立网站存放文件夹，一般为hugo的文件夹中<br />执行`hugo new site yourblog`，此文件夹为网站根目录

## 2. 下载主题
hugo主题网站：[https://themes.gohugo.io/](https://themes.gohugo.io/)<br />进入网站根目录中的·themes·文件夹，执行`git clone 主题的GitHub地址`
> 去除主题.git，如果使用过子模块方式，如何删除
> `git rm --cached submodulesname`
> `rm .gitmodules`

## 3. 主题适应性修改
将themes/ananke/exampleSite内的内容，复制到网站根目录下，替换，并修改相应文件内的内容，如网站标题，作者信息等。

## 4. 第一篇博客文章
在网站根目录运行`hugo new post/first.md`。然后使用编辑器编写内容，注意文件头不要更改。完成文章写作后，可使用`hugo server -D`启动本地网站预览结果，`-D`为预览包含草稿（draft：true）。正式发布时，将markdown文件头的draft：true改为false。

## 5. 编译网站
执行`hugo`即可，<br />


# 部署网站
**注意：在部署之前，将网站根目录下的配置文件config.toml中的`baseURL = "**[**https://xinli0111.gitee.io/**](https://xinli0111.gitee.io/)**"`中的网址更新为部署网站网址，否则会出现主题没法显示的情况。**

## github + netlify
> push后netlify自动部署

首先，创建GitHub pages，然后将本地代码push到GitHub。<br />然后，使用GitHub账户登录netlify，创建新网站，选择github仓库。选择部署应用为hugo，部署文件夹为public。

## gitee
> push后需要手动部署下

在服务中，开启page服务，手动部署public文件夹。<br />

