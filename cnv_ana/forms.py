from flask_wtf import FlaskForm
from flask_wtf.file import FileAllowed, FileRequired
from wtforms import Form
from wtforms import SubmitField, StringField, IntegerField, SelectField, SelectMultipleField
from wtforms import widgets
from wtforms.validators import DataRequired, Length, Optional

class MultiCheckboxField(SelectMultipleField):
    widget = widgets.TableWidget(with_table_tag=False)
    option_widget = widgets.CheckboxInput()

    # uncomment to see how the process call passes through this object
    # def process_data(self, value):
    #     return super(MultiCheckboxField, self).process_data(value)

class BasesLookupForm(FlaskForm):
    submit = SubmitField(u'数据分析')
    search = SubmitField(u'结果查询')
    bases = MultiCheckboxField(u'选择基线：')
    batch = SelectField(u'数据批次：')
    outdir = StringField(u'输出目录：', validators=[Optional()])
    parallel = IntegerField(u'并行数量：', validators=[Optional()])
    thread = IntegerField(u'分析线程：', validators=[Optional()])

class AnaForm(FlaskForm):
    submit_ana = SubmitField(u'开始分析')


class ReportForm(FlaskForm):
    submit = SubmitField(u'报告生成')
    search = SubmitField(u'报告查询')
    batch = SelectField(u'数据批次：')
    outdir = StringField(u'输出目录：', validators=[Optional()])

class getReportForm(FlaskForm):
    submit_ana = SubmitField(u'生成报告')